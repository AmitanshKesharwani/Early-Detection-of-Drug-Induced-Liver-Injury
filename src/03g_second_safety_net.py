"""
03g_second_safety_net.py
========================
Script for Step 3, 4, 5 outputs:
3. Print inscope_names.txt in chunks of 250 names, separated by " ; " in the printed output.
4. Second safety net for in-scope rows:
   (a) mol_weight between 700 and 1000
   (b) rows with 3 or more amide bonds on the largest fragment, using RDKit SMARTS [CX3](=O)[NX3]
   Columns: drug_name, title, formula, mol_weight, label, amide_count. Separate lists.
5. In-scope totals: all, labeled, positives, negatives, positive rate, ambiguous.
"""

from pathlib import Path
import pandas as pd
from rdkit import Chem

PROCESSED_DIR = Path("data/processed")
RESOLVED_PATH = PROCESSED_DIR / "structures_resolved.csv"
INSCOPE_PATH = PROCESSED_DIR / "inscope_names.txt"

df_res = pd.read_csv(RESOLVED_PATH)
in_scope_df = df_res[df_res["in_scope"] == True].copy()

# ==============================================================================
# 3. Print inscope_names.txt in chunks of 250 with " ; "
# ==============================================================================
print("=" * 80)
print("TASK 3: IN-SCOPE DRUG NAMES (chunks of 250 separated by ' ; ')")
print("=" * 80)

names = INSCOPE_PATH.read_text(encoding="utf-8").splitlines()
names = [n.strip() for n in names if n.strip()]
print(f"Total in-scope names in file: {len(names)}\n")

chunk_size = 250
num_chunks = (len(names) + chunk_size - 1) // chunk_size
for i in range(0, len(names), chunk_size):
    chunk = names[i : i + chunk_size]
    chunk_idx = i // chunk_size + 1
    print(f"--- Chunk {chunk_idx}/{num_chunks} (names {i+1} to {min(i+chunk_size, len(names))}) ---")
    print(" ; ".join(chunk))
    print()

# ==============================================================================
# 4. Second safety net
# ==============================================================================
print("=" * 80)
print("TASK 4: SECOND SAFETY NET")
print("=" * 80)

amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")

def compute_amide_count(smiles: str) -> int:
    if not isinstance(smiles, str) or not smiles.strip():
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        return 0
    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    matches = largest.GetSubstructMatches(amide_pattern)
    return len(matches)

in_scope_df["amide_count"] = in_scope_df["smiles"].apply(compute_amide_count)

cols = ["drug_name", "title", "formula", "mol_weight", "label", "amide_count"]

# (a) mol_weight between 700 and 1000
cond_a = in_scope_df[(in_scope_df["mol_weight"] >= 700) & (in_scope_df["mol_weight"] <= 1000)].copy()
cond_a = cond_a.sort_values("mol_weight", ascending=False).reset_index(drop=True)

print(f"\n--- 4(a) IN-SCOPE ROWS WITH MOL_WEIGHT BETWEEN 700 AND 1000 ({len(cond_a)} rows) ---")
# Format nicely
print(f"{'drug_name':<32} {'title':<32} {'formula':<22} {'mol_weight':>10} {'label':>6} {'amides':>6}")
print("-" * 112)
for _, r in cond_a.iterrows():
    lbl_str = f"{int(r['label'])}" if pd.notna(r['label']) else "ambig"
    title_str = str(r['title'])[:31] if pd.notna(r['title']) else ""
    formula_str = str(r['formula']) if pd.notna(r['formula']) else ""
    mw_str = f"{r['mol_weight']:.2f}" if pd.notna(r['mol_weight']) else ""
    print(f"{str(r['drug_name']):<32} {title_str:<32} {formula_str:<22} {mw_str:>10} {lbl_str:>6} {r['amide_count']:>6}")

# (b) rows with 3 or more amide bonds on largest fragment
cond_b = in_scope_df[in_scope_df["amide_count"] >= 3].copy()
cond_b = cond_b.sort_values(["amide_count", "mol_weight"], ascending=[False, False]).reset_index(drop=True)

print(f"\n--- 4(b) IN-SCOPE ROWS WITH >= 3 AMIDE BONDS ON LARGEST FRAGMENT ({len(cond_b)} rows) ---")
print(f"{'drug_name':<32} {'title':<32} {'formula':<22} {'mol_weight':>10} {'label':>6} {'amides':>6}")
print("-" * 112)
for _, r in cond_b.iterrows():
    lbl_str = f"{int(r['label'])}" if pd.notna(r['label']) else "ambig"
    title_str = str(r['title'])[:31] if pd.notna(r['title']) else ""
    formula_str = str(r['formula']) if pd.notna(r['formula']) else ""
    mw_str = f"{r['mol_weight']:.2f}" if pd.notna(r['mol_weight']) else ""
    print(f"{str(r['drug_name']):<32} {title_str:<32} {formula_str:<22} {mw_str:>10} {lbl_str:>6} {r['amide_count']:>6}")


# ==============================================================================
# 5. In-scope totals
# ==============================================================================
print("\n" + "=" * 80)
print("TASK 5: UPDATED IN-SCOPE TOTALS")
print("=" * 80)

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

