"""Parse and generate labels from DILIrank 2.0 raw dataset.

Gate checks and rules:
- Verifies raw file integrity via SHA-256.
- Parses sheet 'version 2' skipping title banner.
- Validates and normalizes 'vDILI-Concern' classes.
- Asserts expected category counts (vmost: 217, vless: 351, vno: 414, ambiguous: 354).
- Checks for duplicate CompoundNames and LTKBIDs.
- Outputs labeled dataset (568 positives, 414 negatives = 982 rows) to data/processed/labels.csv.
- Outputs ambiguous drugs (354 rows) to data/processed/ambiguous_drugs.csv.
- Displays class balance and cross-tab of Comment by label.
"""

import hashlib
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile
import pandas as pd

EXPECTED_SHA256 = "1CA1352FF727AF68E68E250EAE2ED775BCA8492335140AC0AFD2233248694993"


def compute_sha256(file_path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(8192):
            hasher.update(chunk)
    return hasher.hexdigest()


def parse_xlsx_sheet(file_path: Path, sheet_name: str = "version 2", skip_rows: int = 1) -> pd.DataFrame:
    """Parse an Excel sheet directly using standard library zipfile and xml parser."""
    with zipfile.ZipFile(file_path) as z:
        # Load shared strings table if present
        shared_strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            sstree = ET.fromstring(z.read("xl/sharedStrings.xml"))
            ns = {"ns": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            for si in sstree.findall("ns:si", ns):
                text = "".join([t.text or "" for t in si.findall(".//ns:t", ns)])
                shared_strings.append(text)

        # Locate sheet file from workbook.xml and rels
        wbtree = ET.fromstring(z.read("xl/workbook.xml"))
        ns = {"ns": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        r_ns = {"r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
        sheets = [
            (s.attrib["name"], s.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"])
            for s in wbtree.findall(".//ns:sheet", ns)
        ]
        sheet_map = dict(sheets)
        if sheet_name not in sheet_map:
            raise ValueError(f"Sheet '{sheet_name}' not found in workbook. Available: {list(sheet_map.keys())}")

        rels_tree = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rns = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
        rel_map = {r.attrib["Id"]: r.attrib["Target"] for r in rels_tree.findall("r:Relationship", rns)}

        target = rel_map[sheet_map[sheet_name]]
        if not target.startswith("xl/"):
            target = "xl/" + target.lstrip("/")

        stree = ET.fromstring(z.read(target))

        rows_by_idx = {}
        max_col = 0

        def col2num(col_str: str) -> int:
            num = 0
            for c in col_str:
                num = num * 26 + (ord(c) - ord("A")) + 1
            return num - 1

        for row in stree.findall(".//ns:row", ns):
            r_num = int(row.attrib.get("r", 0))
            row_dict = {}
            for c in row.findall("ns:c", ns):
                ref = c.attrib.get("r", "")
                m = re.match(r"([A-Z]+)(\d+)", ref)
                if not m:
                    continue
                col_idx = col2num(m.group(1))
                if col_idx > max_col:
                    max_col = col_idx
                t = c.attrib.get("t")
                v = c.find("ns:v", ns)
                val = v.text if v is not None else ""
                if t == "s" and val != "":
                    val = shared_strings[int(val)]
                elif t == "inlineStr":
                    is_elem = c.find(".//ns:t", ns)
                    val = is_elem.text if is_elem is not None else ""
                row_dict[col_idx] = val
            rows_by_idx[r_num] = row_dict

        sorted_rnums = sorted(rows_by_idx.keys())
        grid = [[rows_by_idx[r_num].get(i, "") for i in range(max_col + 1)] for r_num in sorted_rnums]

        headers = grid[skip_rows]
        data = grid[skip_rows + 1 :]

        non_empty_indices = [i for i, h in enumerate(headers) if h.strip()]
        trimmed_headers = [headers[i].strip() for i in non_empty_indices]
        trimmed_data = [[row[i] for i in non_empty_indices] for row in data]

        return pd.DataFrame(trimmed_data, columns=trimmed_headers)


def main():
    raw_dir = Path("data/raw")
    processed_dir = Path("data/processed")
    processed_dir.mkdir(parents=True, exist_ok=True)

    file_name = "Drug Induced Liver Injury Rank (DILIrank 2.0) Dataset  FDA.xlsx"
    file_path = raw_dir / file_name

    print(f"--- 1. Verifying SHA-256 for: {file_name} ---")
    if not file_path.exists():
        raise FileNotFoundError(f"Raw data file not found: {file_path}")

    actual_hash = compute_sha256(file_path)
    print(f"Actual SHA-256:   {actual_hash}")
    print(f"Expected SHA-256: {EXPECTED_SHA256}")
    assert actual_hash.lower() == EXPECTED_SHA256.lower(), (
        f"SHA-256 mismatch! Got {actual_hash}, expected {EXPECTED_SHA256}"
    )
    print("[PASS] SHA-256 hash verified successfully.\n")

    print("--- 2. Loading Sheet 'version 2' ---")
    df = parse_xlsx_sheet(file_path, sheet_name="version 2", skip_rows=1)
    print(f"Loaded {len(df)} rows with columns: {list(df.columns)}")

    print("\n--- 3. Pre-cleaning Raw Value Distributions ---")
    print("\nRaw 'vDILI-Concern' unique values & counts:")
    raw_concern_counts = df["vDILI-Concern"].value_counts(dropna=False)
    print(raw_concern_counts.to_string())

    print("\nRaw 'Comment' unique values & counts:")
    raw_comment_counts = df["Comment"].value_counts(dropna=False)
    print(raw_comment_counts.to_string())

    print("\n--- 4. Cleaning & Normalization ---")
    # 1. Strip whitespace from all text columns
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip()

    # 2. Normalize vDILI-Concern by lowercasing
    df["original_class_raw"] = df["vDILI-Concern"]
    df["original_class"] = df["vDILI-Concern"].str.lower()

    allowed_classes = {
        "vmost-dili-concern",
        "vless-dili-concern",
        "vno-dili-concern",
        "ambiguous-dili-concern",
    }
    actual_classes = set(df["original_class"].unique())
    print(f"Normalized classes found: {actual_classes}")
    assert actual_classes == allowed_classes, (
        f"Unexpected normalized classes: {actual_classes.difference(allowed_classes)}"
    )
    print("[PASS] Normalized classes match expected set exactly.")

    # 3. Assert counts
    concern_counts = df["original_class"].value_counts().to_dict()
    print(f"Normalized class counts: {concern_counts}")

    assert concern_counts.get("vmost-dili-concern", 0) == 217, (
        f"Expected 217 vmost, got {concern_counts.get('vmost-dili-concern')}"
    )
    assert concern_counts.get("vless-dili-concern", 0) == 351, (
        f"Expected 351 vless, got {concern_counts.get('vless-dili-concern')}"
    )
    assert concern_counts.get("vno-dili-concern", 0) == 414, (
        f"Expected 414 vno, got {concern_counts.get('vno-dili-concern')}"
    )
    assert concern_counts.get("ambiguous-dili-concern", 0) == 354, (
        f"Expected 354 ambiguous, got {concern_counts.get('ambiguous-dili-concern')}"
    )
    assert len(df) == 1336, f"Expected 1336 total rows, got {len(df)}"
    print("[PASS] Category counts verified (217 vmost, 351 vless, 414 vno, 354 ambiguous = 1336 total).")

    # 4. Check duplicates in CompoundName (case-insensitive, whitespace-stripped) and in LTKBID
    compound_lower = df["CompoundName"].str.lower()
    dup_compounds = df[compound_lower.duplicated(keep=False)]
    if not dup_compounds.empty:
        print(f"\n[ERROR] Duplicate CompoundNames detected:\n{dup_compounds[['LTKBID', 'CompoundName']]}")
        raise ValueError(f"Duplicate CompoundNames found: {len(dup_compounds)} rows affected.")
    print("Zero duplicate CompoundNames found.")

    dup_ltkbids = df[df["LTKBID"].duplicated(keep=False)]
    if not dup_ltkbids.empty:
        print(f"\n[ERROR] Duplicate LTKBIDs detected:\n{dup_ltkbids[['LTKBID', 'CompoundName']]}")
        raise ValueError(f"Duplicate LTKBIDs found: {len(dup_ltkbids)} rows affected.")
    print("Zero duplicate LTKBIDs found.")
    print("[PASS] No duplicates in CompoundName or LTKBID.")

    print("\n--- 5. Generating Labels & Output Files ---")
    # Standardize column naming
    df = df.rename(
        columns={
            "LTKBID": "ltkb_id",
            "CompoundName": "drug_name",
            "SeverityClass": "severity_class",
            "LabelSection": "label_section",
            "Comment": "comment",
        }
    )

    # 1. Add boolean flag is_new_v2 = (comment == "New")
    df["is_new_v2"] = df["comment"] == "New"

    # Separate ambiguous drugs
    is_ambiguous = df["original_class"] == "ambiguous-dili-concern"
    df_ambiguous = df[is_ambiguous].copy()
    df_labeled = df[~is_ambiguous].copy()

    # Assign binary labels: vmost & vless -> 1, vno -> 0
    df_labeled["label"] = df_labeled["original_class"].map(
        {
            "vmost-dili-concern": 1,
            "vless-dili-concern": 1,
            "vno-dili-concern": 0,
        }
    )
    df_ambiguous["label"] = ""

    # Reorder columns
    target_columns = [
        "ltkb_id",
        "drug_name",
        "label",
        "is_new_v2",
        "original_class_raw",
        "original_class",
        "severity_class",
        "label_section",
        "comment",
    ]

    df_labeled = df_labeled[target_columns]
    df_ambiguous = df_ambiguous[target_columns]

    # Save to data/processed
    labels_csv = processed_dir / "labels.csv"
    ambiguous_csv = processed_dir / "ambiguous_drugs.csv"

    df_labeled.to_csv(labels_csv, index=False)
    df_ambiguous.to_csv(ambiguous_csv, index=False)

    print(f"Saved labeled data to: {labels_csv}")
    print(f"Saved ambiguous data to: {ambiguous_csv}")

    # Assert final dataset counts
    positives = (df_labeled["label"] == 1).sum()
    negatives = (df_labeled["label"] == 0).sum()
    total_labeled = len(df_labeled)
    total_ambiguous = len(df_ambiguous)

    assert positives == 568, f"Expected 568 positives, got {positives}"
    assert negatives == 414, f"Expected 414 negatives, got {negatives}"
    assert total_labeled == 982, f"Expected 982 labeled rows, got {total_labeled}"
    assert total_ambiguous == 354, f"Expected 354 ambiguous rows, got {total_ambiguous}"

    # Assert temporal holdout / new vs older drug counts
    labeled_new = df_labeled[df_labeled["is_new_v2"]]
    labeled_older = df_labeled[~df_labeled["is_new_v2"]]

    assert len(labeled_new) == 169, f"Expected 169 labeled New drugs, got {len(labeled_new)}"
    assert (labeled_new["label"] == 0).sum() == 120, (
        f"Expected 120 New negatives (0), got {(labeled_new['label'] == 0).sum()}"
    )
    assert (labeled_new["label"] == 1).sum() == 49, (
        f"Expected 49 New positives (1), got {(labeled_new['label'] == 1).sum()}"
    )

    assert len(labeled_older) == 813, f"Expected 813 labeled older drugs, got {len(labeled_older)}"
    assert (labeled_older["label"] == 0).sum() == 294, (
        f"Expected 294 older negatives (0), got {(labeled_older['label'] == 0).sum()}"
    )
    assert (labeled_older["label"] == 1).sum() == 519, (
        f"Expected 519 older positives (1), got {(labeled_older['label'] == 1).sum()}"
    )
    print(
        "[PASS] Output assertions verified (568 positives, 414 negatives, 982 total; 354 ambiguous;\n"
        "       New drugs: 169 total [120 neg, 49 pos]; Older drugs: 813 total [294 neg, 519 pos])."
    )

    print("\n--- 6. Summary Statistics & Cross-Tabulations ---")
    print(f"Total Drugs Evaluated: {len(df)}")
    print(f"  - Labeled Drugs:     {total_labeled} (Positives: {positives}, Negatives: {negatives})")
    print(f"  - Ambiguous Drugs:   {total_ambiguous} (Excluded from training/testing)")
    print(f"Class Balance (Positive Rate): {positives / total_labeled:.2%}")

    print("\nCross-tab of Comment by Label (labeled drugs):")
    ctab = pd.crosstab(df_labeled["comment"], df_labeled["label"], margins=True, margins_name="Total")
    ctab.columns = ["Negative (0)", "Positive (1)", "Total"]
    print(ctab.to_string())


if __name__ == "__main__":
    main()

