"""
03c_resolution_review.py
========================
Investigate resolution topics from data/processed/structures_raw.csv:
1. title_mismatch_full (all 55 rows)
2. multimatch_candidates (all CIDs for the 7 multi-match drugs fetched & cached)
3. salt_stripped_four (the 4 rescued salt-stripped drugs)
4. metals_and_inorganics (inorganic drugs without carbon, and metal complexes vs counterion salts)
5. radiotracers (isotope SMILES or isotope name)
6. target_drugs (verification of key benchmark hepatotoxic / target drugs)
7. heavy_non_peptide (MW > 1000, not flagged likely_biologic and not in not_found)
8. fixable_candidates (testing alternative queries for specific not_found drugs)

Outputs saved as review_<topic>.csv in data/processed/.
Does NOT modify structures_raw.csv or any existing dataset.
"""

import difflib
import json
from pathlib import Path
import re
import time
from typing import Any, Dict, List, Optional
import urllib.parse

import pandas as pd
from rdkit import Chem
import requests

# Paths
PROCESSED_DIR = Path("data/processed")
CACHE_MULTI_DIR = Path("data/cache/pubchem_multi")
CACHE_CANDIDATES_DIR = Path("data/cache/pubchem_candidates")

CACHE_MULTI_DIR.mkdir(parents=True, exist_ok=True)
CACHE_CANDIDATES_DIR.mkdir(parents=True, exist_ok=True)

STRUCTURES_PATH = PROCESSED_DIR / "structures_raw.csv"
NAME_REVIEW_PATH = PROCESSED_DIR / "name_review.csv"

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


# ---------------------------------------------------------------------------
# Topic 1: title_mismatch_full
# ---------------------------------------------------------------------------
def run_title_mismatch_full(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("1. TITLE MISMATCH FULL (All rows with ratio < 0.6 after salt stripping)")
    print("=" * 80)

    salt_words = [
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
    salt_re = re.compile(r"\b\w*(?:" + "|".join(salt_words) + r")\b", re.IGNORECASE)

    def clean_name(s):
        if not isinstance(s, str):
            return ""
        s_clean = salt_re.sub("", s.lower())
        return re.sub(r"[^a-z0-9]+", " ", s_clean).strip()

    mismatches = []
    for _, r in df.iterrows():
        drug_name = r["drug_name"]
        title = r["title"]
        if pd.isna(title) or pd.isna(drug_name):
            continue

        c_drug = clean_name(str(drug_name))
        c_title = clean_name(str(title))
        ratio = difflib.SequenceMatcher(None, c_drug, c_title).ratio()

        if ratio < 0.6:
            mismatches.append({
                "drug_name": drug_name,
                "title": title,
                "cid": r["cid"],
                "formula": r["formula"],
                "mol_weight": r["mol_weight"],
                "label": r["label"],
            })

    out_df = pd.DataFrame(mismatches)
    out_path = PROCESSED_DIR / "review_title_mismatch_full.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(out_df)} rows)\n")
    print(out_df.to_string(index=False))
    return out_df


# ---------------------------------------------------------------------------
# Topic 2: multimatch_candidates
# ---------------------------------------------------------------------------
def run_multimatch_candidates(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("2. MULTI-MATCH CANDIDATES (Fetching details for each candidate CID)")
    print("=" * 80)

    mm_rows = df[df["multi_match"] == True]
    records = []

    for _, r in mm_rows.iterrows():
        drug_name = r["drug_name"]
        curr_cid = str(r["cid"]).strip()
        all_cids_str = str(r["all_cids"])
        cid_list = [c.strip() for c in all_cids_str.split(";") if c.strip()]

        # Cap at 30 if more
        if len(cid_list) > 30:
            cid_list = cid_list[:30]

        for cid in cid_list:
            cache_file = CACHE_MULTI_DIR / f"{cid}.json"
            if cache_file.exists():
                try:
                    data = json.loads(cache_file.read_text(encoding="utf-8"))
                except Exception:
                    data = {}
            else:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularFormula,MolecularWeight/JSON"
                resp = get_with_backoff(url)
                if resp.status_code == 200:
                    try:
                        data = resp.json()
                        cache_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
                    except Exception:
                        data = {}
                else:
                    data = {"error": f"HTTP {resp.status_code}"}

            title = None
            formula = None
            mw = None

            if "PropertyTable" in data:
                props = data["PropertyTable"].get("Properties", [])
                if props:
                    p = props[0]
                    title = p.get("Title")
                    formula = p.get("MolecularFormula")
                    mw = p.get("MolecularWeight")

            records.append({
                "drug_name": drug_name,
                "cid": cid,
                "is_currently_used": (cid == curr_cid or str(int(float(cid))) == str(int(float(curr_cid)))),
                "title": title,
                "formula": formula,
                "mol_weight": mw,
            })

    out_df = pd.DataFrame(records)
    out_path = PROCESSED_DIR / "review_multimatch_candidates.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(out_df)} candidate records)\n")
    print(out_df.to_string(index=False))
    return out_df


# ---------------------------------------------------------------------------
# Topic 3: salt_stripped_four
# ---------------------------------------------------------------------------
def run_salt_stripped_four(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("3. SALT-STRIPPED FOUR DRUGS")
    print("=" * 80)

    ss_df = df[df["match_method"] == "salt_stripped"][
        ["drug_name", "query_used", "cid", "title", "formula", "mol_weight", "smiles"]
    ].copy()

    out_path = PROCESSED_DIR / "review_salt_stripped_four.csv"
    ss_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(ss_df)} rows)\n")
    print(ss_df.to_string(index=False))
    return ss_df


# ---------------------------------------------------------------------------
# Topic 4: metals_and_inorganics
# ---------------------------------------------------------------------------
def run_metals_and_inorganics(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("4. METALS AND INORGANICS (True inorganic & metal complexes vs salt counterions)")
    print("=" * 80)

    resolved = df[df["smiles"].notna()].copy()

    METAL_SYMBOLS = {
        "Pt", "Au", "Fe", "Gd", "Cu", "Co", "Zn", "As", "Hg", "Sb", "Bi",
        "Li", "Mn", "Ga", "Tc", "Lu", "Ra", "Rb", "Yb", "Ag", "Al", "Ca",
        "Mg", "K", "Na", "Se", "Zr", "Ba", "Cr", "Ni", "Ti", "V", "Sn", "Pb"
    }
    SALT_COUNTERION_METALS = {"Na", "K", "Ca", "Mg", "Li"}

    metals_and_inorg = []
    metal_counterions = []

    for _, r in resolved.iterrows():
        smiles = r["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        atoms = mol.GetAtoms()
        has_carbon = any(a.GetAtomicNum() == 6 for a in atoms)
        present_elements = set(a.GetSymbol() for a in atoms)
        metals_present = present_elements & METAL_SYMBOLS

        if not has_carbon:
            metals_and_inorg.append({
                "drug_name": r["drug_name"],
                "label": r["label"],
                "formula": r["formula"],
                "smiles": r["smiles"],
                "reason": "(b) no carbon (inorganic)",
            })
        elif metals_present:
            coord_metals = metals_present - SALT_COUNTERION_METALS
            if coord_metals:
                m_str = ", ".join(sorted(coord_metals))
                metals_and_inorg.append({
                    "drug_name": r["drug_name"],
                    "label": r["label"],
                    "formula": r["formula"],
                    "smiles": r["smiles"],
                    "reason": f"(a) metal/metalloid ({m_str})",
                })
            else:
                m_str = ", ".join(sorted(metals_present))
                metal_counterions.append({
                    "drug_name": r["drug_name"],
                    "label": r["label"],
                    "formula": r["formula"],
                    "smiles": r["smiles"],
                    "reason": f"salt counterion ({m_str})",
                })

    df_inorg = pd.DataFrame(metals_and_inorg)
    out_path = PROCESSED_DIR / "review_metals_and_inorganics.csv"
    df_inorg.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(df_inorg)} rows)\n")
    print("--- TABLE 1: Inorganics (no carbon) & Metal/Metalloid Complexes ---")
    print(df_inorg[["drug_name", "label", "reason", "formula"]].to_string(index=False))

    df_counter = pd.DataFrame(metal_counterions)
    out_path_counter = PROCESSED_DIR / "review_metal_counterions.csv"
    df_counter.to_csv(out_path_counter, index=False)
    print(f"\nSaved: {out_path_counter.resolve()} ({len(df_counter)} rows)")
    print(f"--- TABLE 2: Metal Counterions (Salts in Carbon-containing drugs) --- ({len(df_counter)} drugs)")
    print(df_counter[["drug_name", "label", "reason", "formula"]].head(25).to_string(index=False))
    print(f"... and {len(df_counter) - 25} more salt counterion rows.")
    return df_inorg


# ---------------------------------------------------------------------------
# Topic 5: radiotracers
# ---------------------------------------------------------------------------
def run_radiotracers(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("5. RADIOTRACERS (Isotope bracket in SMILES or isotope token in name)")
    print("=" * 80)

    # Isotope bracket in SMILES: e.g. [18F], [11C], [2H], etc.
    iso_smiles = df["smiles"].str.contains(r"\[\d+[A-Za-z]", regex=True, na=False)

    # Isotope token in name: e.g. f-18, c-11, ga-68, tc-99m, i-123, ra-223, rb-82, yb-169, lu 177
    iso_name = df["drug_name"].str.contains(
        r"(?:f-18|c-11|ga-68|cu-64|tc-99m|i-123|ra-223|rb-82|yb-169|lu\s*177)",
        regex=True,
        case=False,
    )

    rt_df = df[iso_smiles | iso_name][["drug_name", "label", "source", "cid", "smiles"]].copy()
    out_path = PROCESSED_DIR / "review_radiotracers.csv"
    rt_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(rt_df)} rows)\n")
    print(rt_df[["drug_name", "label", "source", "cid"]].to_string(index=False))
    return rt_df


# ---------------------------------------------------------------------------
# Topic 6: target_drugs
# ---------------------------------------------------------------------------
def run_target_drugs(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("6. TARGET BENCHMARK DRUGS VERIFICATION")
    print("=" * 80)

    targets = [
        "isoniazid", "rifampin", "rifampicin", "pyrazinamide",
        "methotrexate", "valproic acid", "valproate",
        "carbamazepine", "phenytoin", "ethambutol"
    ]

    records = []
    for t in targets:
        matches = df[df["drug_name"].str.lower().str.contains(t.lower())]
        if len(matches) == 0:
            records.append({
                "target_query": t,
                "drug_name": "NONE FOUND",
                "source": "—",
                "label": "—",
                "match_method": "—",
                "cid": "—",
                "title": "—",
                "mol_weight": "—",
            })
        else:
            for _, r in matches.iterrows():
                records.append({
                    "target_query": t,
                    "drug_name": r["drug_name"],
                    "source": r["source"],
                    "label": r["label"],
                    "match_method": r["match_method"],
                    "cid": r["cid"],
                    "title": r["title"],
                    "mol_weight": r["mol_weight"],
                })

    out_df = pd.DataFrame(records)
    out_path = PROCESSED_DIR / "review_target_drugs.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(out_df)} rows)\n")
    print(out_df.to_string(index=False))
    return out_df


# ---------------------------------------------------------------------------
# Topic 7: heavy_non_peptide
# ---------------------------------------------------------------------------
def run_heavy_non_peptide(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("7. HEAVY NON-PEPTIDE DRUGS (MW > 1000, not flagged likely_biologic, not not_found)")
    print("=" * 80)

    mw = pd.to_numeric(df["mol_weight"], errors="coerce")
    heavy = df[
        (mw > 1000)
        & (df["likely_biologic"] == False)
        & (df["match_method"] != "not_found")
    ].copy()

    heavy["smiles_len"] = heavy["smiles"].str.len()
    out_df = heavy[["drug_name", "label", "mol_weight", "smiles_len", "title", "cid"]].copy()

    out_path = PROCESSED_DIR / "review_heavy_non_peptide.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(out_df)} rows)\n")
    print(out_df[["drug_name", "label", "mol_weight", "smiles_len"]].to_string(index=False))
    return out_df


# ---------------------------------------------------------------------------
# Topic 8: fixable_candidates
# ---------------------------------------------------------------------------
def run_fixable_candidates(df: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "=" * 80)
    print("8. FIXABLE CANDIDATES LOOKUP (Alternative queries for not_found names)")
    print("=" * 80)

    candidates = [
        ("Ivermectin", "ivermectin"),
        ("Ivermectin", "22,23-dihydroavermectin B1a"),
        ("Polythiazide", "polythiazide"),
        ("Pafuraidine", "pafuramidine"),
        ("Verteporfin", "verteporfin"),
        ("Tetracycline phosphate complex", "tetracycline"),
        ("Phentermine resin complex", "phentermine"),
    ]

    results = []
    for orig_name, alt_query in candidates:
        safe_q = re.sub(r"[^a-zA-Z0-9]+", "_", alt_query.lower()).strip("_")
        cache_file = CACHE_CANDIDATES_DIR / f"{safe_q}.json"

        if cache_file.exists():
            try:
                data = json.loads(cache_file.read_text(encoding="utf-8"))
            except Exception:
                data = {}
        else:
            encoded = urllib.parse.quote(alt_query, safe="")
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/Title,MolecularFormula,MolecularWeight/JSON"
            resp = get_with_backoff(url)
            if resp.status_code == 200:
                try:
                    data = resp.json()
                    cache_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
                except Exception:
                    data = {}
            else:
                data = {"Fault": {"Code": f"HTTP_{resp.status_code}", "Message": "Not Found"}}
                cache_file.write_text(json.dumps(data, indent=2), encoding="utf-8")

        cid = None
        title = None
        formula = None
        mw = None
        status = "not_found"

        if "PropertyTable" in data:
            props = data["PropertyTable"].get("Properties", [])
            if props:
                p = props[0]
                cid = p.get("CID")
                title = p.get("Title")
                formula = p.get("MolecularFormula")
                mw = p.get("MolecularWeight")
                status = "found"

        results.append({
            "drug_name": orig_name,
            "alt_query": alt_query,
            "status": status,
            "cid": cid,
            "title": title,
            "formula": formula,
            "mol_weight": mw,
        })

    out_df = pd.DataFrame(results)
    out_path = PROCESSED_DIR / "review_fixable_candidates.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path.resolve()} ({len(out_df)} candidate queries)\n")
    print(out_df.to_string(index=False))

    print("\n--- Parent Drug Existence Check in structures_raw.csv ---")
    for parent in ["tetracycline", "phentermine"]:
        exact_matches = df[df["drug_name"].str.lower() == parent]
        sub_matches = df[df["drug_name"].str.lower().str.contains(parent)]
        print(f"\nQuerying parent: '{parent}':")
        print(f"  Exact match count: {len(exact_matches)}")
        print(f"  Substring matches ({len(sub_matches)}):")
        print(sub_matches[["ltkb_id", "drug_name", "source", "label", "match_method", "cid"]].to_string(index=False))

    return out_df


def main():
    df = pd.read_csv(STRUCTURES_PATH)
    nr = pd.read_csv(NAME_REVIEW_PATH)

    if "likely_biologic" not in df.columns:
        df = df.merge(nr[["ltkb_id", "likely_biologic"]], on="ltkb_id", how="left")

    run_title_mismatch_full(df)
    run_multimatch_candidates(df)
    run_salt_stripped_four(df)
    run_metals_and_inorganics(df)
    run_radiotracers(df)
    run_target_drugs(df)
    run_heavy_non_peptide(df)
    run_fixable_candidates(df)


if __name__ == "__main__":
    main()

