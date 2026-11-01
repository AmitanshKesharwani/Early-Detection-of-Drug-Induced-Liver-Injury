"""
03d_resolve_and_scope.py
========================
Executes resolution and scoping pipeline:
Part A: Load data/processed/aliases.csv.
Part B: Build data/processed/structures_resolved.csv with alias overrides.
Part C: Generate data/processed/scope_exclusions.csv with precedence rules and add in_scope & exclusion_category columns.
Part E: Print detailed resolution & scoping audit reports.
"""

import json
from pathlib import Path
import re
import time
from typing import Any, Dict, List, Optional
import urllib.parse

import pandas as pd
import requests

PROCESSED_DIR = Path("data/processed")
CACHE_ALIAS_DIR = Path("data/cache/pubchem_alias")
CACHE_ALIAS_DIR.mkdir(parents=True, exist_ok=True)

STRUCTURES_RAW_PATH = PROCESSED_DIR / "structures_raw.csv"
STRUCTURES_RESOLVED_PATH = PROCESSED_DIR / "structures_resolved.csv"
ALIASES_PATH = PROCESSED_DIR / "aliases.csv"
SCOPE_EXCLUSIONS_PATH = PROCESSED_DIR / "scope_exclusions.csv"
RADIOTRACERS_PATH = PROCESSED_DIR / "review_radiotracers.csv"
METALS_PATH = PROCESSED_DIR / "review_metals_and_inorganics.csv"

# Rate limiting
MIN_INTERVAL = 1.0 / 3.0
_last_request_time: float = 0.0


def get_with_backoff(url: str, max_retries: int = 3) -> requests.Response:
    global _last_request_time
    now = time.time()
    elapsed = now - _last_request_time
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)

    delay = 1.0
    for attempt in range(max_retries + 1):
        _last_request_time = time.time()
        try:
            resp = requests.get(url, timeout=15)
        except requests.RequestException:
            if attempt < max_retries:
                time.sleep(delay)
                delay *= 2
                continue
            raise

        if resp.status_code in (429, 500, 502, 503, 504) and attempt < max_retries:
            time.sleep(delay)
            delay *= 2
            _last_request_time = time.time()
            continue

        return resp
    return resp


def fetch_cid_properties(cid: int) -> Dict[str, Any]:
    cache_file = CACHE_ALIAS_DIR / f"{cid}.json"
    if cache_file.exists():
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            pass

    props = "CanonicalSMILES,ConnectivitySMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,Title"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{props}/JSON"
    resp = get_with_backoff(url)
    if resp.status_code == 200:
        data = resp.json()
        cache_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
        return data
    return {}


def main():
    print("=" * 80)
    print("STARTING RESOLUTION AND SCOPING PIPELINE")
    print("=" * 80)

    # -----------------------------------------------------------------------
    # Part A: Read aliases.csv
    # -----------------------------------------------------------------------
    print("\n--- PART A: Loading aliases.csv ---")
    df_aliases = pd.read_csv(ALIASES_PATH)
    print(f"Loaded {len(df_aliases)} rows from {ALIASES_PATH.resolve()}")

    # -----------------------------------------------------------------------
    # Part B: Build structures_resolved.csv
    # -----------------------------------------------------------------------
    print("\n--- PART B: Building structures_resolved.csv ---")
    df_raw = pd.read_csv(STRUCTURES_RAW_PATH)
    df_res = df_raw.copy()

    df_res["resolution"] = "raw"
    df_res["resolution_note"] = ""

    # Build alias lookup map
    alias_map = {}
    for _, a_row in df_aliases.iterrows():
        d_name = str(a_row["drug_name"]).strip()
        alias_map[d_name] = a_row

    changed_rows = []

    for idx, row in df_res.iterrows():
        name = row["drug_name"]
        if name in alias_map:
            a_info = alias_map[name]
            new_cid = int(a_info["cid_override"])
            a_note = str(a_info["note"])

            old_cid = row["cid"]
            old_title = row["title"]
            old_formula = row["formula"]
            old_mw = row["mol_weight"]

            # Fetch properties for new_cid
            prop_data = fetch_cid_properties(new_cid)
            props_list = prop_data.get("PropertyTable", {}).get("Properties", [])
            if props_list:
                p = props_list[0]
                new_title = p.get("Title")
                new_formula = p.get("MolecularFormula")
                new_mw_raw = p.get("MolecularWeight")
                new_mw = float(new_mw_raw) if new_mw_raw is not None else None
                new_smiles = p.get("SMILES") or p.get("IsomericSMILES") or p.get("CanonicalSMILES")
                new_conn_smiles = p.get("ConnectivitySMILES")
                new_inchikey = p.get("InChIKey")

                df_res.at[idx, "cid"] = new_cid
                df_res.at[idx, "title"] = new_title
                df_res.at[idx, "formula"] = new_formula
                df_res.at[idx, "mol_weight"] = new_mw
                df_res.at[idx, "smiles"] = new_smiles
                df_res.at[idx, "connectivity_smiles"] = new_conn_smiles
                df_res.at[idx, "inchikey"] = new_inchikey
                df_res.at[idx, "http_status"] = 200
                df_res.at[idx, "error"] = ""
                df_res.at[idx, "match_method"] = "alias"
                df_res.at[idx, "resolution"] = "alias"
                df_res.at[idx, "resolution_note"] = a_note
                df_res.at[idx, "needs_manual_review"] = False

                changed_rows.append({
                    "drug_name": name,
                    "old_cid": old_cid,
                    "old_title": old_title,
                    "old_formula": old_formula,
                    "old_mw": old_mw,
                    "new_cid": new_cid,
                    "new_title": new_title,
                    "new_formula": new_formula,
                    "new_mw": new_mw,
                })

    print(f"Applied alias overrides to {len(changed_rows)} drugs.")

    # -----------------------------------------------------------------------
    # Part C: Scope Exclusions
    # -----------------------------------------------------------------------
    print("\n--- PART C: Programmatic Scope Exclusions ---")
    radiotracers_df = pd.read_csv(RADIOTRACERS_PATH)
    metals_df = pd.read_csv(METALS_PATH)

    # 1. Radiotracers except Deutetrabenazine
    rad_names = set(radiotracers_df[radiotracers_df["drug_name"].str.lower() != "deutetrabenazine"]["drug_name"])

    # 2. Metals / Inorganics table 1 + Sterile water, Radium ra-223 dichloride, Edetate calcium disodium, Sodium bicarbonate, Lanthanum carbonate
    metal_names = set(metals_df["drug_name"]) | {
        "Sterile water", "Radium ra-223 dichloride", "Edetate calcium disodium",
        "Sodium bicarbonate", "Lanthanum carbonate"
    }

    # 3. Mixture with no single structure (Verteporfin, Estrogens, conjugated, Perflutren, Tetrabamate)
    mixture_names = {"Verteporfin", "Estrogens, conjugated", "Perflutren", "Tetrabamate"}

    # 4. Macromolecules, biologics, polymers
    macro_explicit = {
        "Erythropoietin", "Heparin sodium", "Polidocanol", "Pentosan polysulfate sodium",
        "Aprotinin", "Lepirudin", "Dasiglucagon hydrochloride", "Dalbavancin hydrochloride",
        "Epoetin alfa", "Interferon alfa-2b", "Colistimethate sodium", "Urofollitropin",
        "Eptifibatide", "Pentagastrin"
    }
    # Still not_found in df_res (excluding Verteporfin which has its own mixture category)
    still_not_found = set(df_res[df_res["smiles"].isna()]["drug_name"]) - mixture_names
    macro_names = macro_explicit | still_not_found

    exclusions = []
    exclusion_category_map = {}

    for _, r in df_res.iterrows():
        ltkb_id = r["ltkb_id"]
        name = r["drug_name"]

        reasons = []
        custom_notes = []

        # Precedence:
        # 1. diagnostic_radiotracer
        if name in rad_names:
            reasons.append("diagnostic_radiotracer")

        # 2. inorganic_or_metal_complex
        if name in metal_names:
            reasons.append("inorganic_or_metal_complex")
            if name in {"Sodium bicarbonate", "Lanthanum carbonate"}:
                custom_notes.append("inorganic salt; no C-H or C-C bonds")

        # 3. mixture_no_single_structure
        if name in mixture_names:
            reasons.append("mixture_no_single_structure")
            if name == "Estrogens, conjugated":
                custom_notes.append("mixture; PubChem returned a single component")
            elif name == "Perflutren":
                custom_notes.append("microsphere formulation; C3F8 is not the drug as administered")
            elif name == "Tetrabamate":
                custom_notes.append("mixture of phenobarbital, febarbamate and difebarbamate")

        # 4. macromolecule_biologic_polymer
        if name in macro_names:
            reasons.append("macromolecule_biologic_polymer")
            if name in {"Epoetin alfa", "Interferon alfa-2b"}:
                custom_notes.append("PubChem returned a non-representative model record")
            elif name in {"Colistimethate sodium", "Urofollitropin"}:
                custom_notes.append("PubChem returned a partial record; real molecule is a polypeptide or protein")
            elif name in {"Eptifibatide", "Pentagastrin"}:
                custom_notes.append("peptide of 5 or more residues below 1000 Da; the weight rule cannot catch it")

        if reasons:
            primary = reasons[0]
            other_notes = custom_notes + [f"also {x}" for x in reasons[1:]]
            note_str = "; ".join(other_notes) if other_notes else ""
            exclusions.append({
                "ltkb_id": ltkb_id,
                "drug_name": name,
                "source": r["source"],
                "reason_category": primary,
                "decided_on": "2026-11-01",
                "note": note_str,
                "label": r["label"],
            })
            exclusion_category_map[ltkb_id] = primary

    df_exclusions = pd.DataFrame(exclusions)

    # Assert no duplicate ltkb_id
    assert df_exclusions["ltkb_id"].duplicated().sum() == 0, "Duplicate ltkb_id in scope_exclusions.csv!"
    print(f"Assert passed: no duplicate ltkb_id in {len(df_exclusions)} exclusions.")

    # Save scope_exclusions.csv without internal label column
    cols_to_save = ["ltkb_id", "drug_name", "source", "reason_category", "decided_on", "note"]
    df_exclusions[cols_to_save].to_csv(SCOPE_EXCLUSIONS_PATH, index=False)
    print(f"Written: {SCOPE_EXCLUSIONS_PATH.resolve()} ({len(df_exclusions)} rows)")

    # Add in_scope and exclusion_category to structures_resolved.csv
    df_res["in_scope"] = ~df_res["ltkb_id"].isin(exclusion_category_map)
    df_res["exclusion_category"] = df_res["ltkb_id"].map(exclusion_category_map).fillna("")

    # Add boolean column diagnostic_agent
    diagnostic_list = [
        "Brilliant blue g", "Indocyanine green", "Fluorescein sodium", "Diatrizoate sodium",
        "Iothalamate meglumine", "Ioversol", "Iodipamide sodium", "Aminohippurate sodium",
        "Histamine phosphate", "Regadenoson", "Macimorelin acetate", "Protirelin",
        "Edrophonium chloride", "Sincalide", "Pentagastrin", "Cosyntropin"
    ]
    diag_lower_set = {d.lower() for d in diagnostic_list}
    df_res["diagnostic_agent"] = df_res["drug_name"].str.lower().isin(diag_lower_set)

    # Check diagnostic agent matches
    res_names_lower = set(df_res["drug_name"].str.lower())
    matched_diag = [d for d in diagnostic_list if d.lower() in res_names_lower]
    not_found_diag = [d for d in diagnostic_list if d.lower() not in res_names_lower]
    print(f"\nDiagnostic agents matched ({len(matched_diag)} of {len(diagnostic_list)}): {matched_diag}")
    if not_found_diag:
        print(f"Diagnostic agents NOT found ({len(not_found_diag)}): {not_found_diag}")
    else:
        print("All diagnostic agents matched.")

    # Cast cid to Int64
    if "cid" in df_res.columns:
        df_res["cid"] = df_res["cid"].astype("Int64")

    df_res.to_csv(STRUCTURES_RESOLVED_PATH, index=False)
    print(f"Written: {STRUCTURES_RESOLVED_PATH.resolve()} ({len(df_res)} rows)")

    # -----------------------------------------------------------------------
    # Part E: Print Audit & Summary Outputs
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("E.1: ROWS CHANGED BY ALIAS")
    print("=" * 80)
    df_ch = pd.DataFrame(changed_rows)
    for _, cr in df_ch.iterrows():
        print(f"\nDrug: {cr['drug_name']}")
        print(f"  OLD: CID={cr['old_cid']}, Title={repr(cr['old_title'])}, Formula={cr['old_formula']}, MW={cr['old_mw']}")
        print(f"  NEW: CID={cr['new_cid']}, Title={repr(cr['new_title'])}, Formula={cr['new_formula']}, MW={cr['new_mw']}")

    print("\n" + "=" * 80)
    print("E.2: TABLE OF EXCLUSIONS BY CATEGORY AND LABEL")
    print("=" * 80)
    crosstab = pd.crosstab(
        df_exclusions["reason_category"],
        df_exclusions["label"].fillna("ambiguous"),
        margins=True,
    )
    print(crosstab.to_string())

    print("\n" + "=" * 80)
    print("E.3: IN-SCOPE COUNTS COMPUTED FROM DATA")
    print("=" * 80)
    in_scope_df = df_res[df_res["in_scope"] == True]
    total_in_scope = len(in_scope_df)
    labeled_in_scope = in_scope_df[in_scope_df["source"] == "labeled"]
    ambiguous_in_scope = in_scope_df[in_scope_df["source"] == "ambiguous"]
    pos_count = (labeled_in_scope["label"] == 1).sum()
    neg_count = (labeled_in_scope["label"] == 0).sum()
    pos_rate = (pos_count / len(labeled_in_scope)) * 100

    print(f"Total in-scope drugs     : {total_in_scope}")
    print(f"Labeled in-scope drugs   : {len(labeled_in_scope)}")
    print(f"  - Positive (label 1)   : {pos_count}")
    print(f"  - Negative (label 0)   : {neg_count}")
    print(f"  - In-scope positive rate: {pos_rate:.2f}%")
    print(f"Ambiguous in-scope drugs : {len(ambiguous_in_scope)}")

    print("\n" + "=" * 80)
    print("E.4: TARGET BENCHMARK DRUGS IN-SCOPE CHECK")
    print("=" * 80)
    target_list = [
        "Isoniazid",
        "Rifampin",
        "Pyrazinamide",
        "Methotrexate sodium",
        "Valproic acid",
        "Carbamazepine",
        "Phenytoin",
        "Ethambutol hydrochloride",
    ]
    for target in target_list:
        match = df_res[df_res["drug_name"] == target]
        if len(match) == 0:
            print(f"CRITICAL ERROR: Target drug '{target}' NOT FOUND in structures_resolved.csv!")
            raise SystemExit(1)
        in_s = bool(match.iloc[0]["in_scope"])
        print(f"  {target:<26}: in_scope = {in_s}")
        if not in_s:
            print(f"CRITICAL STOP: Target drug '{target}' is excluded!")
            raise SystemExit(1)

    print("\nAll target benchmark drugs verified as in_scope == True.")

    print("\n" + "=" * 80)
    print("E.5: EXPANDED SAFETY-NET LIST")
    print("=" * 80)
    expanded_terms = [
        "insulin", "globulin", "toxin", "alfa", "beta-", "gamma", "recombinant",
        "peg", "hormone", "factor", "vaccine", "antigen", "extract", "serum",
        "interferon", "epoetin", "poetin", "grastim", "stim", "ermin", "plase",
        "tropin", "albumin", "hyaluron", "heparin", "glucagon", "calcitonin", "oxin"
    ]
    pattern = "|".join(expanded_terms)

    # Search in drug_name or title
    name_match = in_scope_df["drug_name"].str.contains(pattern, case=False, na=False)
    title_match = in_scope_df["title"].str.contains(pattern, case=False, na=False)
    safety_matches = in_scope_df[name_match | title_match]
    print(f"Matches found: {len(safety_matches)}\n")
    cols_to_print = ["drug_name", "title", "formula", "mol_weight", "label"]
    print(safety_matches[cols_to_print].to_string(index=False))

    # -----------------------------------------------------------------------
    # Part 8: Write inscope_names.txt
    # -----------------------------------------------------------------------
    in_scope_names = sorted(in_scope_df["drug_name"].dropna().unique().tolist())
    inscope_file = PROCESSED_DIR / "inscope_names.txt"
    inscope_file.write_text("\n".join(in_scope_names), encoding="utf-8")
    print(f"\nWritten {len(in_scope_names)} in-scope drug names to {inscope_file.resolve()}")


if __name__ == "__main__":
    main()
