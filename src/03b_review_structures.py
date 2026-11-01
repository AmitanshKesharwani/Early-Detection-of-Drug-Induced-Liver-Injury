"""
03b_review_structures.py
========================
Generate audit and review CSV files from data/processed/structures_raw.csv.

Does NOT modify structures_raw.csv or any raw data file.

Outputs written to data/processed/:
1. review_random20.csv
2. review_multimatch.csv
3. review_title_mismatch.csv
4. review_heavy.csv
5. review_not_found.csv
"""

import difflib
from pathlib import Path
import re
import pandas as pd

# Paths
PROCESSED_DIR = Path("data/processed")
STRUCTURES_PATH = PROCESSED_DIR / "structures_raw.csv"
NAME_REVIEW_PATH = PROCESSED_DIR / "name_review.csv"

OUT_RANDOM20 = PROCESSED_DIR / "review_random20.csv"
OUT_MULTIMATCH = PROCESSED_DIR / "review_multimatch.csv"
OUT_TITLE_MISMATCH = PROCESSED_DIR / "review_title_mismatch.csv"
OUT_HEAVY = PROCESSED_DIR / "review_heavy.csv"
OUT_NOT_FOUND = PROCESSED_DIR / "review_not_found.csv"

# Salt / hydrate words to strip for fuzzy comparison
SALT_WORDS = [
    "sulfate", "hydrochloride", "sodium", "potassium", "calcium",
    "maleate", "mesylate", "tartrate", "acetate", "phosphate",
    "citrate", "besylate", "bromide", "chloride", "succinate",
    "fumarate", "hydrobromide", "hydrate", "malate", "gluconate",
    "propionate", "nitrate", "iodide", "valerate", "tosylate",
    "carbonate", "bicarbonate", "lactate", "oxalate", "benzoate",
    "salicylate", "decanoate", "enanthate", "cypionate", "pivalate",
    "tromethamine", "meglumine", "disodium", "trisodium", "dipotassium",
    "monohydrate", "dihydrate", "trihydrate", "dihydrochloride"
]
SALT_REGEX = re.compile(r"\b\w*(?:" + "|".join(SALT_WORDS) + r")\b", re.IGNORECASE)


def clean_name_for_compare(text: str) -> str:
    """Lowercase, strip salt/hydrate terms and non-alphanumeric chars."""
    if not isinstance(text, str):
        return ""
    cleaned = SALT_REGEX.sub("", text.lower())
    cleaned = re.sub(r"[^a-z0-9]+", " ", cleaned).strip()
    return cleaned


def main():
    print("=" * 80)
    print("RUNNING STRUCTURES AUDIT & REVIEW")
    print("=" * 80)

    # 1. Load data
    df = pd.read_csv(STRUCTURES_PATH)
    nr = pd.read_csv(NAME_REVIEW_PATH)

    # Merge likely_biologic if needed
    if "likely_biologic" not in df.columns:
        df = df.merge(nr[["ltkb_id", "likely_biologic"]], on="ltkb_id", how="left")

    # Format cid as integer where present
    if "cid" in df.columns:
        df["cid"] = df["cid"].astype("Int64")

    # -----------------------------------------------------------------------
    # 1. review_random20.csv: 20 random rows with match_method == exact (seed 42)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("1. REVIEW RANDOM 20 (match_method == 'exact', random_state=42)")
    print("=" * 80)
    exact_df = df[df["match_method"] == "exact"]
    r20 = exact_df.sample(n=20, random_state=42)[["drug_name", "title", "cid", "mol_weight"]].copy()
    r20.to_csv(OUT_RANDOM20, index=False)
    print(f"Saved: {OUT_RANDOM20.resolve()} ({len(r20)} rows)\n")
    print(r20.to_string(index=False))

    # -----------------------------------------------------------------------
    # 2. review_multimatch.csv: every row with multi_match == True
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("2. REVIEW MULTI-MATCH (multi_match == True)")
    print("=" * 80)
    mm_df = df[df["multi_match"] == True][["drug_name", "n_cids", "all_cids", "cid", "title"]].copy()
    mm_df.rename(columns={"cid": "cid used"}, inplace=True)
    mm_df.to_csv(OUT_MULTIMATCH, index=False)
    print(f"Saved: {OUT_MULTIMATCH.resolve()} ({len(mm_df)} rows)\n")
    print(mm_df.to_string(index=False))

    # -----------------------------------------------------------------------
    # 3. review_title_mismatch.csv: difflib ratio < 0.6
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("3. REVIEW TITLE MISMATCH (ratio < 0.6 after salt stripping)")
    print("=" * 80)
    mismatch_rows = []
    for _, row in df.iterrows():
        drug_name = row["drug_name"]
        title = row["title"]
        if pd.isna(title) or pd.isna(drug_name):
            continue

        c_drug = clean_name_for_compare(str(drug_name))
        c_title = clean_name_for_compare(str(title))
        ratio = difflib.SequenceMatcher(None, c_drug, c_title).ratio()

        if ratio < 0.6:
            mismatch_rows.append({
                "drug_name": drug_name,
                "title": title,
                "cid": row["cid"],
                "match_method": row["match_method"],
            })

    title_mismatch_df = pd.DataFrame(mismatch_rows)
    title_mismatch_df.to_csv(OUT_TITLE_MISMATCH, index=False)
    print(f"Saved: {OUT_TITLE_MISMATCH.resolve()} ({len(title_mismatch_df)} rows)\n")
    print(title_mismatch_df.to_string(index=False))

    # -----------------------------------------------------------------------
    # 4. review_heavy.csv: mol_weight > 1000 and count table by bins
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("4. REVIEW HEAVY (mol_weight > 1000)")
    print("=" * 80)
    mw_numeric = pd.to_numeric(df["mol_weight"], errors="coerce")
    heavy_df = df[mw_numeric > 1000][["drug_name", "source", "label", "mol_weight", "likely_biologic", "match_method"]].copy()
    heavy_df.to_csv(OUT_HEAVY, index=False)
    print(f"Saved: {OUT_HEAVY.resolve()} ({len(heavy_df)} rows)\n")
    print(heavy_df.to_string(index=False))

    # Print MW count table in bins: <200, 200-500, 500-1000, 1000-2000, >2000
    print("\n--- MOLECULAR WEIGHT DISTRIBUTION TABLE ---")
    bins = [-float("inf"), 200, 500, 1000, 2000, float("inf")]
    bin_labels = ["<200", "200-500", "500-1000", "1000-2000", ">2000"]
    binned_series = pd.cut(mw_numeric, bins=bins, labels=bin_labels, right=False)
    bin_counts = binned_series.value_counts()[bin_labels]

    mw_summary = pd.DataFrame({
        "MW Bin": bin_labels,
        "Drug Count": [bin_counts[lbl] for lbl in bin_labels],
        "Percentage": [f"{(bin_counts[lbl] / mw_numeric.notna().sum()) * 100:.1f}%" for lbl in bin_labels],
    })
    print(mw_summary.to_string(index=False))
    print(f"\nTotal with molecular weight: {mw_numeric.notna().sum()}")
    print(f"Not found / missing MW     : {mw_numeric.isna().sum()}")

    # -----------------------------------------------------------------------
    # 5. review_not_found.csv: every not_found row
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("5. REVIEW NOT FOUND (109 drugs)")
    print("=" * 80)
    nf_df = df[df["match_method"] == "not_found"][["ltkb_id", "drug_name", "source", "label", "likely_biologic", "match_method", "http_status", "error"]].copy()
    nf_df.to_csv(OUT_NOT_FOUND, index=False)
    print(f"Saved: {OUT_NOT_FOUND.resolve()} ({len(nf_df)} rows)\n")

    nf_bio_false = nf_df[nf_df["likely_biologic"] == False]
    print(f"--- NOT FOUND with likely_biologic == False ({len(nf_bio_false)} drugs) ---")
    print(nf_bio_false[["drug_name", "source", "label", "error"]].to_string(index=False))

    nf_bio_true = nf_df[nf_df["likely_biologic"] == True]
    print(f"\n--- NOT FOUND with likely_biologic == True ({len(nf_bio_true)} drugs) ---")
    print(", ".join(sorted(nf_bio_true["drug_name"].tolist())))

    # -----------------------------------------------------------------------
    # 6. Fragment ('.') and Isotope bracket counts in SMILES
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("6. SMILES FRAGMENTS AND ISOTOPE COUNTS")
    print("=" * 80)
    smiles_series = df["smiles"].dropna()

    # Multi-fragment: contains '.'
    dot_count = smiles_series.str.contains(r"\.", regex=True).sum()

    # Isotope bracket: contains '[' followed by digits and element symbol, e.g. [18F], [11C], [2H]
    iso_count = smiles_series.str.contains(r"\[\d+[A-Za-z]", regex=True).sum()

    print(f"Total resolved structures (with SMILES): {len(smiles_series)}")
    print(f"Rows with '.' in SMILES (multi-fragment/salts): {dot_count} ({dot_count / len(smiles_series) * 100:.1f}%)")
    print(f"Rows with isotope bracket in SMILES (e.g. [18F]): {iso_count} ({iso_count / len(smiles_series) * 100:.1f}%)")


if __name__ == "__main__":
    main()

