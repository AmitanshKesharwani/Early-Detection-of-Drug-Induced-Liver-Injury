import pandas as pd

raw = pd.read_csv("data/processed/structures_raw.csv")
res = pd.read_csv("data/processed/structures_resolved.csv")
excl = pd.read_csv("data/processed/scope_exclusions.csv")

old_row = raw[raw["drug_name"] == "Dextromethorphan polistirex"].iloc[0]
new_row = res[res["drug_name"] == "Dextromethorphan polistirex"].iloc[0]

print("=" * 75)
print("RECORD CHANGED BY DEXTROMETHORPHAN ALIAS (OLD vs NEW)")
print("=" * 75)
print("Drug Name   :", old_row["drug_name"])
print(f"  OLD RECORD: CID={old_row['cid']}, Title={repr(old_row['title'])}, Formula={old_row['formula']}, MW={old_row['mol_weight']}")
print(f"  NEW RECORD: CID={new_row['cid']}, Title={repr(new_row['title'])}, Formula={new_row['formula']}, MW={new_row['mol_weight']}")
print("  SMILES    :", new_row["smiles"])
print(f"  Resolution: {new_row['resolution']}, Note: {new_row['resolution_note']}")

print("\n" + "=" * 75)
print("EXCLUSIONS BY CATEGORY AND LABEL")
print("=" * 75)
excl_with_label = excl.merge(res[["ltkb_id", "label"]], on="ltkb_id", how="left")
ct_detail = pd.crosstab(
    excl_with_label["reason_category"],
    excl_with_label["label"].map({0.0: "neg (0)", 1.0: "pos (1)"}).fillna("ambiguous"),
    margins=True,
)
print(ct_detail.to_string())

print("\n" + "=" * 75)
print("IN-SCOPE TOTALS COMPUTED FROM THE DATA")
print("=" * 75)
in_scope = res[res["in_scope"] == True]
labeled_in_scope = in_scope[in_scope["source"] == "labeled"]
ambig_in_scope = in_scope[in_scope["source"] == "ambiguous"]
pos = (labeled_in_scope["label"] == 1).sum()
neg = (labeled_in_scope["label"] == 0).sum()
pos_rate = (pos / len(labeled_in_scope)) * 100

print(f"  Total in-scope drugs     : {len(in_scope)}")
print(f"  Labeled in-scope drugs   : {len(labeled_in_scope)}")
print(f"    - Positives (label 1)  : {pos}")
print(f"    - Negatives (label 0)  : {neg}")
print(f"    - Positive rate        : {pos_rate:.2f}%")
print(f"  Ambiguous in-scope drugs : {len(ambig_in_scope)}")
print(f"  Total excluded drugs     : {len(excl)}")

