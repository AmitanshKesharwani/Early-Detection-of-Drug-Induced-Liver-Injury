"""
02b_pubchem_pilot.py
====================
Pilot PubChem PUG REST lookup for 9 representative drug names.

Purpose
-------
Demonstrate the PubChem API before writing the full-download script:
  - rate limiting (<=3 req/sec)
  - retry with exponential backoff on 429/5xx
  - disk caching of every raw response
  - not-found (404/Fault) handling
  - multi-CID awareness
  - special cases: salt SMILES ("." separator), isotope labels, biologics, URL encoding

This is a read-only pilot. Nothing is written to data/processed.

Usage
-----
    python src/02b_pubchem_pilot.py

Output
------
  - Console: per-name details + summary table
  - data/cache/pubchem_pilot/<safe_name>_cids.json    (step-1 CID list)
  - data/cache/pubchem_pilot/<safe_name>_props.json   (step-2 properties)
"""

import json
import pathlib
import re
import time
import urllib.parse
import requests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CACHE_DIR = pathlib.Path("data/cache/pubchem_pilot")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"

# Properties we request (PubChem silently maps old names to new ones).
PROPS = "CanonicalSMILES,ConnectivitySMILES,IsomericSMILES,InChIKey,MolecularFormula,MolecularWeight"

# Rate limiting: <=3 requests per second.
MIN_INTERVAL = 1.0 / 3.0   # ~0.333 s between requests

# ---------------------------------------------------------------------------
# Pilot lookup list
# Each dict: query_name (sent to API), display_name (shown in output), note
# ---------------------------------------------------------------------------
LOOKUPS = [
    {"query": "aspirin",               "display": "aspirin",               "note": ""},
    {"query": "isoniazid",             "display": "isoniazid",             "note": ""},
    {"query": "methotrexate",          "display": "methotrexate",          "note": ""},
    {"query": "acetaminophen",         "display": "acetaminophen",         "note": ""},
    {"query": "valproic acid",         "display": "valproic acid",         "note": ""},
    {"query": "Abacavir sulfate",      "display": "Abacavir sulfate",      "note": "CHECK: SMILES contains '.'?"},
    {"query": "abacavir",              "display": "abacavir (bare)",       "note": "EXTRA: compare with sulfate form"},
    {"query": "Abatacept",             "display": "Abatacept",             "note": "EXPECTED 404 (biologic)"},
    {"query": "Florbetaben f-18",      "display": "Florbetaben f-18",      "note": "Isotope-labeled name"},
    {"query": "florbetaben",           "display": "florbetaben (bare)",    "note": "EXTRA: compare with isotope form"},
    {"query": "Estrogens, conjugated", "display": "Estrogens, conjugated", "note": "Combination name; URL-encode comma"},
]

# ---------------------------------------------------------------------------
# Helper: safe filename for cache
# ---------------------------------------------------------------------------

def safe_filename(name: str) -> str:
    """Convert a query name to a filesystem-safe string."""
    s = name.lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


# ---------------------------------------------------------------------------
# Helper: rate-limited GET with retry
# ---------------------------------------------------------------------------

_last_request_time: float = 0.0


def get_with_backoff(url: str, max_retries: int = 3) -> requests.Response:
    """
    GET url, honouring the <=3 req/sec rate limit.
    Retries on 429 (Too Many Requests) or 5xx with exponential backoff.
    """
    global _last_request_time

    # Rate limiting: sleep if we are going too fast.
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
            print(f"    [Network error, attempt {attempt+1}] {exc}")
            if attempt < max_retries:
                print(f"    Sleeping {delay}s before retry ...")
                time.sleep(delay)
                delay *= 2
                continue
            raise

        if resp.status_code in (429, 500, 502, 503, 504) and attempt < max_retries:
            print(f"    [HTTP {resp.status_code}, attempt {attempt+1}] Sleeping {delay}s before retry ...")
            time.sleep(delay)
            delay *= 2
            # Also honour the rate limit after retry sleep.
            _last_request_time = time.time()
            continue

        return resp

    return resp   # last attempt result (may be error status)


# ---------------------------------------------------------------------------
# Core lookup: one query name -> {cid_info, prop_info}
# ---------------------------------------------------------------------------

def lookup_name(query: str) -> dict:
    """
    Perform two-step PubChem lookup for *query*.

    Step 1: /compound/name/{encoded}/cids/JSON   -> list of all CIDs
    Step 2: /compound/name/{encoded}/property/{PROPS}/JSON -> properties (best match)

    Returns a dict with keys:
        query, encoded, cids_url, props_url,
        cids_status, cids_json,
        props_status, props_json,
        n_cids, cid_used,
        smiles_key, smiles, connectivity_smiles,
        inchikey, mol_formula, mol_weight,
        raw_prop_keys
    """
    encoded = urllib.parse.quote(query, safe="")
    cids_url = f"{BASE_URL}/{encoded}/cids/JSON"
    props_url = f"{BASE_URL}/{encoded}/property/{PROPS}/JSON"

    safe = safe_filename(query)
    cids_cache = CACHE_DIR / f"{safe}_cids.json"
    props_cache = CACHE_DIR / f"{safe}_props.json"

    result = {
        "query": query,
        "encoded": encoded,
        "cids_url": cids_url,
        "props_url": props_url,
    }

    # ---- Step 1: CID list ----
    if cids_cache.exists():
        print(f"    [cache] {cids_cache.name}")
        cids_text = cids_cache.read_text(encoding="utf-8")
        cids_json = json.loads(cids_text)
        cids_status = 200  # cached = was successful
    else:
        resp = get_with_backoff(cids_url)
        cids_status = resp.status_code
        cids_text = resp.text
        cids_json = None
        try:
            cids_json = resp.json()
        except Exception:
            cids_json = {"_raw": cids_text}
        cids_cache.write_text(json.dumps(cids_json, indent=2), encoding="utf-8")

    result["cids_status"] = cids_status
    result["cids_json"] = cids_json

    # Parse CID list.
    n_cids = 0
    all_cids = []
    if isinstance(cids_json, dict):
        if "IdentifierList" in cids_json:
            all_cids = cids_json["IdentifierList"].get("CID", [])
            n_cids = len(all_cids)
        # Fault response -> n_cids stays 0

    result["n_cids"] = n_cids
    result["all_cids"] = all_cids

    # ---- Step 2: Properties (best-match) ----
    if props_cache.exists():
        print(f"    [cache] {props_cache.name}")
        props_text = props_cache.read_text(encoding="utf-8")
        props_json = json.loads(props_text)
        props_status = 200
    else:
        resp2 = get_with_backoff(props_url)
        props_status = resp2.status_code
        props_text = resp2.text
        props_json = None
        try:
            props_json = resp2.json()
        except Exception:
            props_json = {"_raw": props_text}
        props_cache.write_text(json.dumps(props_json, indent=2), encoding="utf-8")

    result["props_status"] = props_status
    result["props_json"] = props_json

    # Parse properties.
    raw_prop_keys = []
    cid_used = None
    smiles_key = None
    smiles = None
    conn_smiles = None
    inchikey = None
    mol_formula = None
    mol_weight = None

    if isinstance(props_json, dict) and "PropertyTable" in props_json:
        props_list = props_json["PropertyTable"].get("Properties", [])
        if props_list:
            p = props_list[0]
            raw_prop_keys = list(p.keys())
            cid_used = p.get("CID")
            # SMILES key is literally "SMILES" in current PubChem responses
            if "SMILES" in p:
                smiles_key = "SMILES"
                smiles = p["SMILES"]
            elif "IsomericSMILES" in p:
                smiles_key = "IsomericSMILES"
                smiles = p["IsomericSMILES"]
            elif "CanonicalSMILES" in p:
                smiles_key = "CanonicalSMILES"
                smiles = p["CanonicalSMILES"]
            conn_smiles = p.get("ConnectivitySMILES")
            inchikey = p.get("InChIKey")
            mol_formula = p.get("MolecularFormula")
            mol_weight = p.get("MolecularWeight")

    result["raw_prop_keys"] = raw_prop_keys
    result["cid_used"] = cid_used
    result["smiles_key"] = smiles_key
    result["smiles"] = smiles
    result["connectivity_smiles"] = conn_smiles
    result["inchikey"] = inchikey
    result["mol_formula"] = mol_formula
    result["mol_weight"] = mol_weight

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    sep = "=" * 72
    thin = "-" * 60

    print(sep)
    print("PubChem PUG REST pilot -- 9 representative names")
    print(f"Cache directory: {CACHE_DIR.resolve()}")
    print(sep)

    results_by_display = {}

    for entry in LOOKUPS:
        display = entry["display"]
        query   = entry["query"]
        note    = entry["note"]

        print(f"\n{thin}")
        print(f"NAME    : {display}")
        if note:
            print(f"NOTE    : {note}")
        print(f"QUERY   : {repr(query)}")

        r = lookup_name(query)
        results_by_display[display] = r

        # ---- Step 1 output ----
        print(f"\n[Step 1 -- CID list]")
        print(f"  URL    : {r['cids_url']}")
        print(f"  Status : {r['cids_status']}")
        print(f"  n_CIDs : {r['n_cids']}")
        if r["n_cids"] > 1:
            print(f"  All CIDs (first 5): {r['all_cids'][:5]}")
        if r["n_cids"] == 0:
            # Show fault structure
            print(f"  Fault JSON: {json.dumps(r['cids_json'], indent=4)}")

        # ---- Step 2 output ----
        print(f"\n[Step 2 -- Properties]")
        print(f"  URL    : {r['props_url']}")
        print(f"  Status : {r['props_status']}")
        if r["raw_prop_keys"]:
            print(f"  Raw JSON keys returned: {r['raw_prop_keys']}")
        else:
            print(f"  Fault/empty JSON: {json.dumps(r['props_json'], indent=4)}")

        if r["cid_used"] is not None:
            print(f"  CID used          : {r['cid_used']}")
        if r["smiles_key"]:
            print(f"  SMILES key        : {r['smiles_key']}")
            smiles_preview = (r["smiles"] or "")[:120]
            print(f"  SMILES (preview)  : {smiles_preview}")
            has_dot = "." in (r["smiles"] or "")
            print(f"  SMILES has '.'?   : {has_dot}  {'<-- multi-fragment/salt' if has_dot else ''}")
        if r["connectivity_smiles"]:
            cs_preview = (r["connectivity_smiles"] or "")[:80]
            print(f"  ConnectivitySMILES: {cs_preview}")
        if r["inchikey"]:
            print(f"  InChIKey          : {r['inchikey']}")
        if r["mol_formula"]:
            print(f"  MolecularFormula  : {r['mol_formula']}")
        if r["mol_weight"]:
            print(f"  MolecularWeight   : {r['mol_weight']}")

    # -----------------------------------------------------------------------
    # Special investigative cases
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("SPECIAL INVESTIGATION RESULTS")
    print(sep)

    # --- Florbetaben f-18 vs florbetaben (bare) ---
    r_f18  = results_by_display.get("Florbetaben f-18")
    r_bare = results_by_display.get("florbetaben (bare)")
    print("\n[Florbetaben: isotope vs bare]")
    if r_f18 and r_bare:
        cid_f18  = r_f18["cid_used"]
        cid_bare = r_bare["cid_used"]
        same_cid = (cid_f18 == cid_bare)
        print(f"  Florbetaben f-18  CID : {cid_f18}")
        print(f"  florbetaben (bare) CID: {cid_bare}")
        print(f"  Same CID?             : {same_cid}")
        smiles_f18 = r_f18.get("smiles") or ""
        has_18f = "[18F]" in smiles_f18
        print(f"  IsomericSMILES of f-18 contains [18F]? : {has_18f}")
        cs_f18 = r_f18.get("connectivity_smiles") or ""
        plain_f = ("F" in cs_f18) and ("[18F]" not in cs_f18)
        print(f"  ConnectivitySMILES uses plain F (no isotope)? : {plain_f}")

    # --- Abacavir sulfate vs abacavir (bare) ---
    r_salt = results_by_display.get("Abacavir sulfate")
    r_free = results_by_display.get("abacavir (bare)")
    print("\n[Abacavir: salt vs free base]")
    if r_salt and r_free:
        cid_salt = r_salt["cid_used"]
        cid_free = r_free["cid_used"]
        same_cid = (cid_salt == cid_free)
        print(f"  Abacavir sulfate CID : {cid_salt}")
        print(f"  abacavir (bare)  CID : {cid_free}")
        print(f"  Same CID?            : {same_cid}")
        smiles_salt = r_salt.get("smiles") or ""
        has_dot_salt = "." in smiles_salt
        smiles_free = r_free.get("smiles") or ""
        has_dot_free = "." in smiles_free
        print(f"  Salt SMILES has '.' (multi-fragment)? : {has_dot_salt}")
        print(f"  Free-base SMILES has '.'?             : {has_dot_free}")
        if has_dot_salt:
            fragments = smiles_salt.split(".")
            print(f"  Salt SMILES fragment count: {len(fragments)}")

    # --- Abatacept (biologic, expected 404) ---
    r_bio = results_by_display.get("Abatacept")
    print("\n[Abatacept -- expected 404 (biologic)]")
    if r_bio:
        print(f"  CID step status  : {r_bio['cids_status']}")
        print(f"  Prop step status : {r_bio['props_status']}")
        print(f"  n_CIDs           : {r_bio['n_cids']}")
        print(f"  Fault JSON       : {json.dumps(r_bio['cids_json'], indent=4)}")

    # --- Estrogens, conjugated (URL-encoding demo) ---
    r_estro = results_by_display.get("Estrogens, conjugated")
    print("\n[Estrogens, conjugated -- URL encoding demo]")
    if r_estro:
        print(f"  Encoded query    : {r_estro['encoded']}")
        print(f"  CID used         : {r_estro['cid_used']}")
        smiles_estro = r_estro.get("smiles") or ""
        print(f"  SMILES has '.'?  : {'.' in smiles_estro}")

    # -----------------------------------------------------------------------
    # NOT-FOUND example (for reference)
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("NOT-FOUND RESPONSE FORMAT (from Abatacept or similar)")
    print(sep)
    for display, r in results_by_display.items():
        cj = r.get("cids_json") or {}
        if isinstance(cj, dict) and "Fault" in cj:
            print(f"\nQuery: {r['query']}")
            print(f"Fault JSON:\n{json.dumps(cj, indent=4)}")
            break

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("SUMMARY TABLE")
    print(sep)
    col_widths = [26, 8, 8, 12, 24, 14]
    headers = ["name", "HTTP", "n_CIDs", "CID used", "SMILES key", "mol weight"]
    header_line = "  ".join(h.ljust(w) for h, w in zip(headers, col_widths))
    print(header_line)
    print("-" * len(header_line))

    for entry in LOOKUPS:
        display = entry["display"]
        r = results_by_display.get(display, {})
        cols = [
            display,
            str(r.get("props_status", "--")),
            str(r.get("n_cids", "--")),
            str(r.get("cid_used", "--")),
            str(r.get("smiles_key") or "--"),
            str(r.get("mol_weight") or "--"),
        ]
        row = "  ".join(c.ljust(w) for c, w in zip(cols, col_widths))
        print(row)

    print(f"\n{sep}")
    print("Cache files written to:", CACHE_DIR.resolve())
    cache_files = sorted(CACHE_DIR.glob("*.json"))
    print(f"  {len(cache_files)} file(s) cached:")
    for f in cache_files:
        print(f"    {f.name}  ({f.stat().st_size} bytes)")
    print(sep)


if __name__ == "__main__":
    main()
