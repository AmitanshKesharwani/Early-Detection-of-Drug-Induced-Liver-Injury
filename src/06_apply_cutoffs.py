"""
Step 6 – Apply the 1000 Da cutoff and collapse duplicates.

Input : data/processed/parents_candidate.csv  (1,167 rows, 1 per in-scope drug)
Output: data/processed/parents_final.csv       (one row per unique drug kept for modelling)

Rules
-----
1. Drop rows where parent_mw > 1000.  These are saved to
   review6_heavy_dropped.csv for the record.

2. Conflict group (different labels on the same connectivity_key):
   Decision (user 2026-11-01): keep BOTH rows with their own labels.
   No rows are removed; both participate in train/test splits as
   separate drugs.

3. Mixed labeled + ambiguous groups (6 groups):
   Decision (user 2026-11-01): keep BOTH rows unchanged.
   The labeled member keeps its label; the ambiguous member stays
   ambiguous (label NaN). No rows are removed or relabelled.

4. Same-label duplicate groups (18 groups, all labeled; or all ambiguous):
   Collapse each group to ONE representative row.
   - Tie-breaking priority (descending):
       a. is_new_v2 == True  (prefer DILIrank 2.0 additions)
       b. source == "DILIrank"  (prefer primary source over DILIst)
       c. shorter drug_name (alphabetical proxy for the "base" salt)
       d. lowest ltkb_id lexicographically (deterministic last resort)
   - All drug_names in the group are recorded in a new column
     `collapsed_from` (pipe-separated) on the representative row.
   - Rows removed by collapse are saved to review6_collapsed.csv.

5. Assertions (gate checks):
   a. No ltkb_id appears twice in parents_final.csv.
   b. No drug appears in both labeled train/test pool and the
      locked test set (not relevant here, but assert no NaN
      connectivity_key).
   c. Label counts before == label counts after for kept labeled rows
      (i.e. no re-labelling happened).
   d. Every row in parents_final.csv had in_scope==True in
      structures_resolved.csv (join on ltkb_id and check).

6. Add a boolean column `use_taut_smiles` = True for all rows
   (records the decision that parent_smiles is parent_smiles_taut).
   Also add `step6_applied` = True so downstream scripts can verify
   they are reading the final file.

After writing, print:
- Row counts in / out / dropped-heavy / collapsed
- Label distribution (labeled total, positives, negatives, ambiguous)
- Positive rate among labeled
- The conflict group (amphetamine) with both rows printed
- The 6 mixed groups (one line each)
- The 18 collapsed groups (old members → representative)
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROCESSED = Path("data/processed")
INPUT = PROCESSED / "parents_candidate.csv"
OUTPUT = PROCESSED / "parents_final.csv"
HEAVY_REVIEW = PROCESSED / "review6_heavy_dropped.csv"
COLLAPSED_REVIEW = PROCESSED / "review6_collapsed.csv"

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT, dtype={"ltkb_id": str})
n_in = len(df)
print(f"=== Step 6: Apply cutoffs and collapse duplicates ===")
print(f"Input rows: {n_in}")
print(f"Label distribution in: pos={int((df['label']==1).sum())}, "
      f"neg={int((df['label']==0).sum())}, ambig={int(df['label'].isna().sum())}")

# ---------------------------------------------------------------------------
# 2. Apply 1000 Da cutoff
# ---------------------------------------------------------------------------
heavy_mask = df["parent_mw"] > 1000
heavy_dropped = df[heavy_mask].copy()
df_after_mw = df[~heavy_mask].copy()

heavy_dropped.to_csv(HEAVY_REVIEW, index=False)
print(f"\n--- 1000 Da cutoff ---")
print(f"Dropped (parent_mw > 1000): {len(heavy_dropped)} rows")
print(f"  pos dropped: {int((heavy_dropped['label']==1).sum())}, "
      f"neg dropped: {int((heavy_dropped['label']==0).sum())}, "
      f"ambig dropped: {int(heavy_dropped['label'].isna().sum())}")
print(f"Remaining after MW cutoff: {len(df_after_mw)}")
for _, row in heavy_dropped.iterrows():
    lbl = "ambig" if pd.isna(row["label"]) else int(row["label"])
    print(f"  HEAVY: {row['drug_name']}  parent_mw={row['parent_mw']:.1f}  label={lbl}")

# ---------------------------------------------------------------------------
# 3. Identify duplicate groups
# ---------------------------------------------------------------------------
dup_keys = (df_after_mw.groupby("connectivity_key")
            .filter(lambda x: len(x) > 1)["connectivity_key"]
            .unique())

print(f"\n--- Duplicate groups after MW cutoff ---")
print(f"Unique connectivity_keys with >=2 rows: {len(dup_keys)}")

# Classify groups
conflict_groups = []     # labeled rows with different labels
mixed_groups = []        # labeled + ambiguous rows (same label for labeled)
same_label_groups = []   # all labeled same label, or all ambiguous

for key in dup_keys:
    grp = df_after_mw[df_after_mw["connectivity_key"] == key]
    labeled_rows = grp[grp["label"].notna()]
    ambig_rows = grp[grp["label"].isna()]
    unique_labels = labeled_rows["label"].unique()

    if len(unique_labels) > 1:
        conflict_groups.append(key)
    elif len(ambig_rows) > 0 and len(labeled_rows) > 0:
        mixed_groups.append(key)
    else:
        same_label_groups.append(key)

print(f"  Conflict groups (different labels): {len(conflict_groups)}")
print(f"  Mixed labeled+ambiguous groups:     {len(mixed_groups)}")
print(f"  Same-label / all-ambiguous groups:  {len(same_label_groups)}")

# ---------------------------------------------------------------------------
# 4. Handle conflict group: keep both rows
# ---------------------------------------------------------------------------
print(f"\n--- Conflict group (keep both rows) ---")
for key in conflict_groups:
    grp = df_after_mw[df_after_mw["connectivity_key"] == key]
    for _, row in grp.iterrows():
        print(f"  CONFLICT key={key}: {row['drug_name']}  label={int(row['label'])}  "
              f"ltkb_id={row['ltkb_id']}")

# ---------------------------------------------------------------------------
# 5. Handle mixed labeled+ambiguous groups: keep both rows
# ---------------------------------------------------------------------------
print(f"\n--- Mixed labeled+ambiguous groups (keep both rows) ---")
for key in mixed_groups:
    grp = df_after_mw[df_after_mw["connectivity_key"] == key]
    names = grp["drug_name"].tolist()
    labels = ["ambig" if pd.isna(l) else int(l) for l in grp["label"]]
    print(f"  MIXED key={key}: {list(zip(names, labels))}")

# ---------------------------------------------------------------------------
# 6. Collapse same-label duplicate groups
# ---------------------------------------------------------------------------
print(f"\n--- Same-label / all-ambiguous groups (collapse) ---")

# Tie-breaking helper: lower = preferred representative
def sort_key(row):
    return (
        0 if row.get("is_new_v2") is True or str(row.get("is_new_v2")).lower() == "true" else 1,
        0 if str(row.get("source", "")).lower() == "dilirank" else 1,
        len(str(row.get("drug_name", ""))),
        str(row.get("ltkb_id", ""))
    )

rows_to_collapse = []  # ltkb_ids of the non-representatives

# Add collapsed_from column (pipe-separated list of all drug_names in group)
df_after_mw = df_after_mw.copy()
df_after_mw["collapsed_from"] = df_after_mw["drug_name"]  # default = self

for key in same_label_groups:
    grp = df_after_mw[df_after_mw["connectivity_key"] == key].copy()
    grp_sorted = grp.apply(sort_key, axis=1)
    grp = grp.iloc[grp_sorted.argsort()]  # sort ascending = best first
    all_names = grp["drug_name"].tolist()
    rep_idx = grp.index[0]
    non_rep_idxs = grp.index[1:].tolist()

    collapsed_from_str = " | ".join(all_names)
    df_after_mw.at[rep_idx, "collapsed_from"] = collapsed_from_str

    rows_to_collapse.extend(non_rep_idxs)

    rep_name = grp.iloc[0]["drug_name"]
    non_rep_names = [grp.loc[i, "drug_name"] for i in non_rep_idxs]
    labels_in = ["ambig" if pd.isna(grp.loc[i, "label"]) else int(grp.loc[i, "label"])
                 for i in grp.index]
    print(f"  COLLAPSE key={key}: {all_names} labels={labels_in} -> rep='{rep_name}'")

# Save collapsed rows for review
collapsed_review_df = df_after_mw.loc[rows_to_collapse].copy()
collapsed_review_df.to_csv(COLLAPSED_REVIEW, index=False)
print(f"\nCollapsed away: {len(rows_to_collapse)} rows -> review6_collapsed.csv")

# Remove non-representative rows from same-label groups
df_final = df_after_mw.drop(index=rows_to_collapse).copy()

# ---------------------------------------------------------------------------
# 7. Add bookkeeping columns
# ---------------------------------------------------------------------------
df_final["use_taut_smiles"] = True
df_final["step6_applied"] = True

# ---------------------------------------------------------------------------
# 8. Assertions / gate checks
# ---------------------------------------------------------------------------
print(f"\n--- Assertions ---")

# a. No duplicate ltkb_id
assert df_final["ltkb_id"].nunique() == len(df_final), \
    "FAIL: duplicate ltkb_id in parents_final.csv"
print("[PASS] No duplicate ltkb_id")

# b. No null connectivity_key
assert df_final["connectivity_key"].notna().all(), \
    "FAIL: null connectivity_key found"
print("[PASS] No null connectivity_key")

# c. No re-labelling: every labeled row has label in {0.0, 1.0}
labeled_final = df_final[df_final["label"].notna()]
assert set(labeled_final["label"].unique()).issubset({0.0, 1.0}), \
    "FAIL: unexpected label values"
print("[PASS] Labels are only 0 or 1 (no relabelling)")

# d. Every row in parents_final was in_scope in structures_resolved
resolved = pd.read_csv(PROCESSED / "structures_resolved.csv", dtype={"ltkb_id": str})
scope_ids = set(resolved[resolved["in_scope"] == True]["ltkb_id"])
final_ids = set(df_final["ltkb_id"])
assert final_ids.issubset(scope_ids), \
    f"FAIL: rows in parents_final not in-scope in structures_resolved: {final_ids - scope_ids}"
print("[PASS] All rows were in_scope in structures_resolved.csv")

# e. Label counts for labeled rows: check none were silently dropped beyond MW cutoff
# Among labeled rows, we only remove: (a) heavy >1000Da, (b) collapsed non-reps with same label
# Any label-count change should only come from (a) + (b)
n_pos_in = int((df["label"] == 1).sum())
n_neg_in = int((df["label"] == 0).sum())
n_ambig_in = int(df["label"].isna().sum())
n_pos_final = int((df_final["label"] == 1).sum())
n_neg_final = int((df_final["label"] == 0).sum())
n_ambig_final = int(df_final["label"].isna().sum())
# Verify heavy drop + collapse explains the difference
heavy_pos = int((heavy_dropped["label"] == 1).sum())
heavy_neg = int((heavy_dropped["label"] == 0).sum())
heavy_ambig = int(heavy_dropped["label"].isna().sum())
collapsed_df = df_after_mw.loc[rows_to_collapse]
coll_pos = int((collapsed_df["label"] == 1).sum())
coll_neg = int((collapsed_df["label"] == 0).sum())
coll_ambig = int(collapsed_df["label"].isna().sum())

assert n_pos_in - heavy_pos - coll_pos == n_pos_final, \
    f"FAIL pos count: {n_pos_in} - {heavy_pos} - {coll_pos} != {n_pos_final}"
assert n_neg_in - heavy_neg - coll_neg == n_neg_final, \
    f"FAIL neg count: {n_neg_in} - {heavy_neg} - {coll_neg} != {n_neg_final}"
assert n_ambig_in - heavy_ambig - coll_ambig == n_ambig_final, \
    f"FAIL ambig count: {n_ambig_in} - {heavy_ambig} - {coll_ambig} != {n_ambig_final}"
print("[PASS] Label counts are consistent (heavy drops + collapses account for all changes)")

# ---------------------------------------------------------------------------
# 9. Write output
# ---------------------------------------------------------------------------
df_final.to_csv(OUTPUT, index=False)
print(f"\n=== parents_final.csv written ===")
print(f"Output rows: {len(df_final)}")
print(f"  Labeled:   {n_pos_final + n_neg_final}")
print(f"    Positive (1): {n_pos_final}")
print(f"    Negative (0): {n_neg_final}")
print(f"  Ambiguous: {n_ambig_final}")
if (n_pos_final + n_neg_final) > 0:
    pos_rate = 100 * n_pos_final / (n_pos_final + n_neg_final)
    print(f"  Positive rate (labeled only): {pos_rate:.1f}%")
print(f"\nDropped by MW cutoff: {len(heavy_dropped)}")
print(f"Collapsed (same-label duplicates removed): {len(rows_to_collapse)}")
print(f"  review6_heavy_dropped.csv: {len(heavy_dropped)} rows")
print(f"  review6_collapsed.csv:     {len(rows_to_collapse)} rows")
