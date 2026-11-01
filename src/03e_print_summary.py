"""
03e_print_summary.py
====================
Print Step 6-8 output:
  6. In-scope totals, exclusions table, full aliases.csv text.
  7. Safety-net expanded (already computed -- re-run from structures_resolved.csv).
  8. inscope_names.txt contents in chunks of 300.
"""
import sys
import textwrap
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
RESOLVED = BASE / "data/processed/structures_resolved.csv"
EXCLUSIONS = BASE / "data/processed/scope_exclusions.csv"
ALIASES = BASE / "data/processed/aliases.csv"
INSCOPE_NAMES = BASE / "data/processed/inscope_names.txt"

# ─── Step 6: in-scope totals ────────────────────────────────────────────────
print("=" * 70)
print("STEP 6 — IN-SCOPE TOTALS AND EXCLUSION TABLE")
print("=" * 70)

res = pd.read_csv(RESOLVED)
excl = pd.read_csv(EXCLUSIONS)
aliases = pd.read_csv(ALIASES)

in_scope = res[res["in_scope"] == True]
out_scope = res[res["in_scope"] == False]

total_inscope = len(in_scope)
labeled_inscope = in_scope[in_scope["source"] == "labeled"]
pos = (labeled_inscope["label"] == 1).sum()
neg = (labeled_inscope["label"] == 0).sum()
ambig_inscope = (in_scope["source"] == "ambiguous").sum()
pos_rate = pos / (pos + neg) * 100 if (pos + neg) > 0 else float("nan")

print(f"\n  Total in-scope       : {total_inscope}")
print(f"  Labeled in-scope     : {len(labeled_inscope)}")
print(f"  Positives (label=1)  : {pos}")
print(f"  Negatives (label=0)  : {neg}")
print(f"  Positive rate        : {pos_rate:.2f}%")
print(f"  Ambiguous in-scope   : {ambig_inscope}")
print(f"  Total excluded       : {len(excl)}")

# Exclusions by category and label
print("\n  Exclusions by category and label:")
excl2 = excl.copy()
excl2 = excl2.merge(res[["ltkb_id", "label"]], on="ltkb_id", how="left")
excl2["label_str"] = excl2["label"].apply(
    lambda x: "pos" if x == 1 else ("neg" if x == 0 else "ambiguous")
)
tbl = (
    excl2.groupby(["reason_category", "label_str"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)
# Ensure columns
for col in ["neg", "pos", "ambiguous"]:
    if col not in tbl.columns:
        tbl[col] = 0
tbl["All"] = tbl["neg"] + tbl["pos"] + tbl["ambiguous"]
tbl = tbl.sort_values("reason_category").reset_index(drop=True)
print(f"\n  {'reason_category':<40} {'neg':>5} {'pos':>5} {'ambig':>7} {'All':>5}")
print("  " + "-" * 62)
for _, row in tbl.iterrows():
    print(f"  {row['reason_category']:<40} {row['neg']:>5} {row['pos']:>5} {row['ambiguous']:>7} {row['All']:>5}")
total_row = tbl[["neg","pos","ambiguous","All"]].sum()
print("  " + "-" * 62)
print(f"  {'All':<40} {int(total_row['neg']):>5} {int(total_row['pos']):>5} {int(total_row['ambiguous']):>7} {int(total_row['All']):>5}")

# Full aliases.csv text
print("\n\n  Full text of aliases.csv:")
print("  " + "-" * 66)
for line in ALIASES.read_text(encoding="utf-8").splitlines():
    print("  " + line)
print("  " + "-" * 66)

# ─── Step 7: expanded safety-net ────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 7 — EXPANDED SAFETY-NET (in-scope rows matching biologic keywords)")
print("=" * 70)

keywords = [
    "insulin", "globulin", "toxin", "alfa", "beta-", "gamma", "recombinant",
    "peg", "hormone", "factor", "vaccine", "antigen", "extract", "serum",
    "interferon", "epoetin", "poetin", "grastim", "stim", "ermin", "plase",
    "tropin", "albumin", "hyaluron", "heparin", "glucagon", "calcitonin",
    "oxin",
]

def contains_any_keyword(text: str) -> bool:
    t = str(text).lower()
    return any(kw.lower() in t for kw in keywords)

in_scope_df = res[res["in_scope"] == True].copy()
# Fill title NaN for matching
in_scope_df["_title"] = in_scope_df["title"].fillna("")
mask = in_scope_df["drug_name"].apply(contains_any_keyword) | in_scope_df["_title"].apply(contains_any_keyword)
hits = in_scope_df[mask][["drug_name", "title", "formula", "mol_weight", "label"]].copy()
hits = hits.sort_values("drug_name").reset_index(drop=True)

print(f"\n  {len(hits)} in-scope rows match safety-net keywords:\n")
if hits.empty:
    print("  (none)")
else:
    col_w = {"drug_name": 30, "title": 30, "formula": 22, "mol_weight": 10, "label": 6}
    header = (
        f"  {'drug_name':<30} {'title':<30} {'formula':<22} {'mol_weight':>10} {'label':>6}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for _, row in hits.iterrows():
        lbl = str(row["label"]).replace(".0", "") if pd.notna(row["label"]) else "ambig"
        print(
            f"  {str(row['drug_name']):<30} "
            f"{str(row['title'])[:29]:<30} "
            f"{str(row['formula']):<22} "
            f"{row['mol_weight']:>10} "
            f"{lbl:>6}"
        )

# ─── Step 8: inscope_names.txt in chunks of 300 ─────────────────────────────
print("\n" + "=" * 70)
print("STEP 8 — inscope_names.txt CONTENTS (chunks of 300)")
print("=" * 70)

names_text = INSCOPE_NAMES.read_text(encoding="utf-8").strip()
names = [n.strip() for n in names_text.split(",") if n.strip()]
print(f"\n  Total names: {len(names)}\n")

CHUNK = 300
for i in range(0, len(names), CHUNK):
    chunk = names[i : i + CHUNK]
    chunk_num = i // CHUNK + 1
    total_chunks = (len(names) + CHUNK - 1) // CHUNK
    print(f"  --- Chunk {chunk_num}/{total_chunks} (names {i+1}-{min(i+CHUNK, len(names))}) ---")
    # Wrap at 100 chars
    line = ", ".join(chunk)
    for wrapped_line in textwrap.wrap(line, width=100):
        print("  " + wrapped_line)
    print()

print("=" * 70)
print("DONE")
print("=" * 70)

