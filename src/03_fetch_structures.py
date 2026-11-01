"""
03_fetch_structures.py
======================
Look up chemical structures on PubChem for all 1,336 drugs in name_review.csv.

Lookup Order (stops at first success):
1. alias: query from data/processed/aliases.csv if present
2. exact: drug name as given
3. isotope_stripped: trailing isotope token (e.g. ' f-18', ' c-11', ' ga-68') removed
4. salt_stripped: trailing salt/hydrate words removed
If all fail: match_method = 'not_found'

API & Cache Rules:
- Rate limit: <= 3 requests per second
- Exponential backoff on 429 and 5xx
- Disk cache at data/cache/pubchem/ (every raw query cached, never queried twice)
- Output: data/processed/structures_raw.csv and data/processed/pubchem_failures.csv
- Supports --limit N for testing (e.g. --limit 25)
"""

import argparse
import hashlib
import json
import pathlib
import re
import time
import urllib.parse
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Paths and Constants
# ---------------------------------------------------------------------------
CACHE_DIR = pathlib.Path("data/cache/pubchem")
NAME_REVIEW_PATH = pathlib.Path("data/processed/name_review.csv")
ALIASES_PATH = pathlib.Path("data/processed/aliases.csv")
OUTPUT_PATH = pathlib.Path("data/processed/structures_raw.csv")
FAILURES_PATH = pathlib.Path("data/processed/pubchem_failures.csv")

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
PROPS = "CanonicalSMILES,ConnectivitySMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight,Title"

MIN_INTERVAL = 1.0 / 3.0  # Max 3 requests per second (~0.333s)
_last_request_time: float = 0.0

# ---------------------------------------------------------------------------
# Stripping patterns
# ---------------------------------------------------------------------------
# Trailing isotope token (e.g., ' f-18', ' c-11', ' ga-68', ' cu-64', ' tc-99m')
ISOTOPE_TRAILING_PATTERN = re.compile(
    r"\s+(?:[a-zA-Z]{1,2}[-\s]\d{1,3}[a-zA-Z]?)$",
    re.IGNORECASE,
)

# Salt / hydrate roots to strip from the end of the name
SALT_ROOTS = [
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
SALT_TRAILING_PATTERN = re.compile(
    r"\s+(?:s-)?\b\w*(?:" + "|".join(SALT_ROOTS) + r")\b$",
    re.IGNORECASE,
)


def strip_isotope(name: str) -> str:
    """Remove trailing isotope token from drug name."""
    stripped = ISOTOPE_TRAILING_PATTERN.sub("", name).strip()
    return stripped if stripped else name


def strip_salt(name: str) -> str:
    """Repeatedly strip trailing salt/hydrate words until stable or single word."""
    curr = name
    while True:
        next_val = SALT_TRAILING_PATTERN.sub("", curr).strip()
        if next_val == curr or not next_val:
            break
        curr = next_val
    return curr if curr else name


# ---------------------------------------------------------------------------
# Cache Helpers
# ---------------------------------------------------------------------------
def get_cache_paths(query: str) -> Tuple[pathlib.Path, pathlib.Path]:
    """Return (cids_cache_path, props_cache_path) for a given query string."""
    safe_slug = re.sub(r"[^a-zA-Z0-9]+", "_", query.lower()).strip("_")[:40]
    q_hash = hashlib.sha256(query.encode("utf-8")).hexdigest()[:10]
    prefix = f"{safe_slug}_{q_hash}"
    return (
        CACHE_DIR / f"{prefix}_cids.json",
        CACHE_DIR / f"{prefix}_props.json",
    )


# ---------------------------------------------------------------------------
# Rate-limited HTTP GET with exponential backoff
# ---------------------------------------------------------------------------
def get_with_backoff(url: str, max_retries: int = 3) -> requests.Response:
    """
    GET url honoring <= 3 req/sec rate limit.
    Retries on 429 and 5xx with exponential backoff.
    """
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
        except requests.RequestException as exc:
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
# Single query lookup against PubChem
# ---------------------------------------------------------------------------
def query_pubchem(query: str) -> Tuple[bool, Dict[str, Any]]:
    """
    Execute 2-step lookup for `query`.
    Returns (success: bool, info: dict).
    """
    encoded = urllib.parse.quote(query, safe="")
    cids_url = f"{BASE_URL}/{encoded}/cids/JSON"
    props_url = f"{BASE_URL}/{encoded}/property/{PROPS}/JSON"

    cids_path, props_path = get_cache_paths(query)

    # Step 1: CIDs
    if cids_path.exists():
        try:
            cids_json = json.loads(cids_path.read_text(encoding="utf-8"))
            cids_status = 200 if "IdentifierList" in cids_json else 404
        except Exception:
            cids_json = {}
            cids_status = 500
    else:
        resp = get_with_backoff(cids_url)
        cids_status = resp.status_code
        try:
            cids_json = resp.json()
        except Exception:
            cids_json = {"_raw": resp.text}
        # Cache deterministic responses (200 OK or 404 Fault)
        if cids_status in (200, 404):
            cids_path.write_text(json.dumps(cids_json, indent=2), encoding="utf-8")

    all_cids = []
    if isinstance(cids_json, dict) and "IdentifierList" in cids_json:
        all_cids = cids_json["IdentifierList"].get("CID", [])

    if not all_cids:
        # Not found or empty CID list
        error_msg = "No CID found"
        if isinstance(cids_json, dict) and "Fault" in cids_json:
            error_msg = cids_json["Fault"].get("Message", "No CID found")
        return False, {
            "query": query,
            "http_status": cids_status,
            "n_cids": 0,
            "cid": None,
            "all_cids": "",
            "multi_match": False,
            "smiles": None,
            "connectivity_smiles": None,
            "inchikey": None,
            "formula": None,
            "mol_weight": None,
            "title": None,
            "error": error_msg,
        }

    # Step 2: Properties
    if props_path.exists():
        try:
            props_json = json.loads(props_path.read_text(encoding="utf-8"))
            props_status = 200
        except Exception:
            props_json = {}
            props_status = 500
    else:
        resp_props = get_with_backoff(props_url)
        props_status = resp_props.status_code
        try:
            props_json = resp_props.json()
        except Exception:
            props_json = {"_raw": resp_props.text}
        if props_status in (200, 404):
            props_path.write_text(json.dumps(props_json, indent=2), encoding="utf-8")

    cid_used = all_cids[0]
    multi_match = len(all_cids) > 1
    all_cids_str = ";".join(str(c) for c in all_cids)

    smiles = None
    connectivity_smiles = None
    inchikey = None
    formula = None
    mol_weight = None
    title = None
    error = None

    if isinstance(props_json, dict) and "PropertyTable" in props_json:
        props_list = props_json["PropertyTable"].get("Properties", [])
        if props_list:
            p = props_list[0]
            cid_used = p.get("CID", cid_used)
            # PubChem property names:
            smiles = p.get("SMILES") or p.get("IsomericSMILES") or p.get("CanonicalSMILES")
            connectivity_smiles = p.get("ConnectivitySMILES")
            inchikey = p.get("InChIKey")
            formula = p.get("MolecularFormula")
            mol_weight = p.get("MolecularWeight")
            title = p.get("Title")
    else:
        error = "Failed to parse property table"

    return True, {
        "query": query,
        "http_status": props_status,
        "n_cids": len(all_cids),
        "cid": cid_used,
        "all_cids": all_cids_str,
        "multi_match": multi_match,
        "smiles": smiles,
        "connectivity_smiles": connectivity_smiles,
        "inchikey": inchikey,
        "formula": formula,
        "mol_weight": mol_weight,
        "title": title,
        "error": error,
    }


# ---------------------------------------------------------------------------
# Drug lookup logic following the 4-step order
# ---------------------------------------------------------------------------
def resolve_drug(drug_name: str, aliases_map: Dict[str, str]) -> Tuple[str, Dict[str, Any]]:
    """
    Try lookup stages:
    1. alias
    2. exact
    3. isotope_stripped
    4. salt_stripped
    Returns (match_method, info_dict).
    """
    # 1. Alias
    if drug_name in aliases_map:
        alias_query = aliases_map[drug_name]
        success, info = query_pubchem(alias_query)
        if success:
            return "alias", info

    # 2. Exact
    success, info = query_pubchem(drug_name)
    if success:
        return "exact", info

    # 3. Isotope stripped
    iso_stripped = strip_isotope(drug_name)
    if iso_stripped != drug_name:
        success, info = query_pubchem(iso_stripped)
        if success:
            return "isotope_stripped", info

    # 4. Salt stripped
    salt_stripped = strip_salt(drug_name)
    if salt_stripped != drug_name:
        success, info = query_pubchem(salt_stripped)
        if success:
            return "salt_stripped", info

    # If all failed: return not_found with info from exact lookup
    return "not_found", info


# ---------------------------------------------------------------------------
# Main Runner
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Fetch chemical structures from PubChem.")
    parser.add_argument("--limit", type=int, default=None, help="Limit to first N drugs.")
    args = parser.parse_args()

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load name_review.csv
    df_review = pd.read_csv(NAME_REVIEW_PATH)
    total_drugs = len(df_review)

    # 2. Load aliases.csv
    aliases_map: Dict[str, str] = {}
    if ALIASES_PATH.exists():
        df_aliases = pd.read_csv(ALIASES_PATH)
        for _, row in df_aliases.iterrows():
            d_name = str(row.get("drug_name", "")).strip()
            q_name = str(row.get("query_name", "")).strip()
            if d_name and q_name and q_name.lower() != "nan":
                aliases_map[d_name] = q_name

    # Apply limit if specified
    if args.limit is not None and args.limit > 0:
        df_target = df_review.head(args.limit).copy()
        print(f"Running in LIMITED mode: first {len(df_target)} of {total_drugs} drugs.")
    else:
        df_target = df_review.copy()
        print(f"Running in FULL mode: all {total_drugs} drugs.")

    print(f"Loaded {len(aliases_map)} alias entries.")
    print("Beginning structure lookups...\n")

    results: List[Dict[str, Any]] = []

    for idx, (_, row) in enumerate(df_target.iterrows()):
        ltkb_id = row["ltkb_id"]
        drug_name = str(row["drug_name"]).strip()
        source = row["source"]
        label_raw = row["label"]
        label = int(label_raw) if pd.notna(label_raw) else ""
        is_new_v2 = bool(row["is_new_v2"])

        method, info = resolve_drug(drug_name, aliases_map)

        query_used = info.get("query", drug_name)
        n_cids = info.get("n_cids", 0)
        cid = info.get("cid")
        all_cids = info.get("all_cids", "")
        multi_match = bool(info.get("multi_match", False))
        smiles = info.get("smiles")
        connectivity_smiles = info.get("connectivity_smiles")
        inchikey = info.get("inchikey")
        formula = info.get("formula")
        mol_weight = info.get("mol_weight")
        title = info.get("title")
        http_status = info.get("http_status")
        error = info.get("error")

        # Evaluate needs_manual_review:
        # Conditions: multi_match, match_method is not exact, combination, likely_biologic, has_digits, mol_weight > 1000, not_found
        comb_flag = bool(row.get("combination", False))
        bio_flag = bool(row.get("likely_biologic", False))
        digits_flag = bool(row.get("has_digits", False))

        high_mw = False
        if mol_weight is not None:
            try:
                if float(mol_weight) > 1000.0:
                    high_mw = True
            except (ValueError, TypeError):
                pass

        needs_review = (
            multi_match
            or (method != "exact")
            or comb_flag
            or bio_flag
            or digits_flag
            or high_mw
            or (method == "not_found")
        )

        results.append({
            "ltkb_id": ltkb_id,
            "drug_name": drug_name,
            "source": source,
            "label": label,
            "is_new_v2": is_new_v2,
            "query_used": query_used,
            "match_method": method,
            "n_cids": n_cids,
            "cid": cid,
            "all_cids": all_cids,
            "multi_match": multi_match,
            "smiles": smiles,
            "connectivity_smiles": connectivity_smiles,
            "inchikey": inchikey,
            "formula": formula,
            "mol_weight": mol_weight,
            "title": title,
            "http_status": http_status,
            "error": error,
            "needs_manual_review": needs_review,
        })

        status_tag = f"CID {cid}" if cid else f"FAIL ({http_status})"
        review_tag = "[REVIEW]" if needs_review else ""
        print(f"[{idx+1:4d}/{len(df_target)}] {drug_name:<30} -> {method:<16} {status_tag:<12} {review_tag}", flush=True)

    df_out = pd.DataFrame(results)
    if "cid" in df_out.columns:
        df_out["cid"] = df_out["cid"].astype("Int64")

    # Save outputs
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved raw structures to: {OUTPUT_PATH.resolve()}")

    df_failures = df_out[df_out["match_method"] == "not_found"].copy()
    df_failures.to_csv(FAILURES_PATH, index=False)
    print(f"Saved {len(df_failures)} failures to: {FAILURES_PATH.resolve()}")

    # Print summary
    print("\n" + "=" * 80)
    print("EXECUTION SUMMARY")
    print("=" * 80)
    print(f"Total processed     : {len(df_out)}")
    print(f"Found on PubChem    : {len(df_out[df_out['match_method'] != 'not_found'])}")
    print(f"Not found (failures): {len(df_failures)}")
    print("\nBreakdown by match_method:")
    for m_method, cnt in df_out["match_method"].value_counts().items():
        print(f"  {m_method:<18}: {cnt:3d}")

    print(f"\nFlagged for manual review: {df_out['needs_manual_review'].sum()} of {len(df_out)}")

    # Print nicely formatted summary table
    print("\n" + "=" * 80)
    print("OUTPUT TABLE PREVIEW (first 25)")
    print("=" * 80)
    cols_to_show = ["drug_name", "match_method", "cid", "mol_weight", "title", "needs_manual_review"]
    print(df_out[cols_to_show].to_string(index=False))
    print("=" * 80)


if __name__ == "__main__":
    main()
