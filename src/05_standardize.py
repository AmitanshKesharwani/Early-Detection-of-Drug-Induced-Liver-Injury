"""
src/05_standardize.py
=====================
Standardizes candidate chemical structures from data/processed/structures_resolved.csv
for all rows with in_scope == True (1,167 rows after Tetrabamate exclusion).

Key Improvements:
1. Heavy-atom-based parent chooser: after Cleanup and Uncharger, splits into fragments,
   discards inorganics (no carbon), and chooses fragment with most heavy atoms.
   Ties flag parent_needs_review = True.
2. Applies parent overrides from data/processed/parent_overrides.csv (Fosfomycin, Methenamine).
3. Explicit hydrogens removed with Chem.RemoveHs after isotope clearance (so deuterated drugs
   like Deutetrabenazine contain no [H]).
4. fragment_flag recomputed on neutralized distinct organic fragments.
5. Tracks parent_old vs parent_new and outputs review5b_parent_changes.csv.
6. Tautomer diagnostic: saves parent_smiles_notaut and parent_smiles_taut, prints tautomer changes
   and imidic acid / oxime formations.
7. Outputs data/processed/parents_candidate.csv (1,167 rows) and review files A-G.
"""

from collections import Counter
import csv
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

PROCESSED_DIR = Path("data/processed")
RESOLVED_PATH = PROCESSED_DIR / "structures_resolved.csv"
OUTPUT_PARENTS_PATH = PROCESSED_DIR / "parents_candidate.csv"
OVERRIDES_PATH = PROCESSED_DIR / "parent_overrides.csv"


def load_overrides() -> Dict[str, Dict[str, str]]:
    """Load parent overrides from data/processed/parent_overrides.csv."""
    if not OVERRIDES_PATH.exists():
        return {}
    df_over = pd.read_csv(OVERRIDES_PATH)
    res = {}
    for _, row in df_over.iterrows():
        res[row["drug_name"]] = {
            "keep_formula": row["keep_formula"],
            "reason": row["reason"],
        }
    return res


def standardize_molecule(
    smiles: str,
    drug_name: str,
    overrides: Dict[str, Dict[str, str]],
    uncharger: rdMolStandardize.Uncharger,
    te: rdMolStandardize.TautomerEnumerator,
    old_lfc: rdMolStandardize.LargestFragmentChooser,
) -> Dict[str, Any]:
    """
    Standardize a single SMILES string using the updated heavy-atom parent selection rule.
    Also computes old chooser parent for comparison.
    """
    res: Dict[str, Any] = {
        "parse_ok": False,
        "n_fragments_before": 0,
        "fragment_flag": False,
        "frag1_smiles": "",
        "frag1_heavy_atoms": 0,
        "frag2_smiles": "",
        "frag2_heavy_atoms": 0,
        "parent_needs_review": False,
        "parent_override_applied": False,
        "parent_smiles_notaut": "",
        "parent_smiles_taut": "",
        "parent_smiles": "",
        "parent_inchikey": "",
        "connectivity_key": "",
        "parent_formula": "",
        "parent_mw": None,
        "heavy_atoms": 0,
        "old_parent_smiles": "",
        "old_connectivity_key": "",
        "old_parent_formula": "",
        "tautomer_changed": False,
        "imidic_acid_formed": False,
    }

    if not isinstance(smiles, str) or not smiles.strip():
        return res

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return res

    res["parse_ok"] = True

    try:
        # ---------------------------------------------------------------------
        # Reference: Compute old chooser result (for comparison)
        # ---------------------------------------------------------------------
        clean_mol_old = rdMolStandardize.Cleanup(mol)
        old_frag = old_lfc.choose(clean_mol_old)
        old_uncharged = uncharger.uncharge(old_frag)
        old_taut = te.Canonicalize(old_uncharged)
        for a in old_taut.GetAtoms():
            if a.GetIsotope() != 0:
                a.SetIsotope(0)
        res["old_parent_smiles"] = Chem.MolToSmiles(old_taut, isomericSmiles=True)
        old_ikey = Chem.MolToInchiKey(old_taut)
        res["old_connectivity_key"] = old_ikey[:14] if old_ikey else ""
        res["old_parent_formula"] = rdMolDescriptors.CalcMolFormula(old_taut)

        # ---------------------------------------------------------------------
        # Step 1: Cleanup and split into fragments
        # ---------------------------------------------------------------------
        clean_mol = rdMolStandardize.Cleanup(mol)
        raw_frags = Chem.GetMolFrags(clean_mol, asMols=True)
        res["n_fragments_before"] = len(raw_frags)

        # Uncharge each fragment and discard inorganic fragments (no carbon)
        org_frags: List[Chem.Mol] = []
        for f in raw_frags:
            if any(a.GetSymbol() == "C" for a in f.GetAtoms()):
                uf = uncharger.uncharge(f)
                org_frags.append(uf)

        if not org_frags:
            org_frags = [uncharger.uncharge(f) for f in raw_frags]

        # ---------------------------------------------------------------------
        # Step 2: Deduplicate distinct fragments after neutralization
        # ---------------------------------------------------------------------
        distinct_dict: Dict[str, Chem.Mol] = {}
        for f in org_frags:
            can_smi = Chem.MolToSmiles(f, isomericSmiles=True)
            if can_smi not in distinct_dict:
                distinct_dict[can_smi] = f

        # Sort distinct fragments descending by heavy-atom count
        sorted_distinct = sorted(
            distinct_dict.values(),
            key=lambda f: f.GetNumHeavyAtoms(),
            reverse=True,
        )

        if len(sorted_distinct) >= 1:
            res["frag1_smiles"] = Chem.MolToSmiles(sorted_distinct[0], isomericSmiles=True)
            res["frag1_heavy_atoms"] = sorted_distinct[0].GetNumHeavyAtoms()

        if len(sorted_distinct) >= 2:
            res["frag2_smiles"] = Chem.MolToSmiles(sorted_distinct[1], isomericSmiles=True)
            res["frag2_heavy_atoms"] = sorted_distinct[1].GetNumHeavyAtoms()
            ha1 = res["frag1_heavy_atoms"]
            ha2 = res["frag2_heavy_atoms"]
            if ha2 > 0 and ha1 <= 1.5 * ha2:
                res["fragment_flag"] = True

        # ---------------------------------------------------------------------
        # Step 3: Choose parent (Overrides -> Heavy-atom count -> Tie detection)
        # ---------------------------------------------------------------------
        chosen_frag: Optional[Chem.Mol] = None

        if drug_name in overrides:
            target_formula = overrides[drug_name]["keep_formula"]
            matches = [
                f for f in sorted_distinct
                if rdMolDescriptors.CalcMolFormula(f) == target_formula
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Override for '{drug_name}' expected exactly 1 fragment with formula "
                    f"'{target_formula}', but found {len(matches)}!"
                )
            chosen_frag = matches[0]
            res["parent_override_applied"] = True
            res["parent_needs_review"] = False
        else:
            if len(sorted_distinct) >= 2:
                ha1 = sorted_distinct[0].GetNumHeavyAtoms()
                ha2 = sorted_distinct[1].GetNumHeavyAtoms()
                if ha1 == ha2:
                    res["parent_needs_review"] = True
            chosen_frag = sorted_distinct[0]

        # ---------------------------------------------------------------------
        # Step 4: Clear isotopes, remove explicit hydrogens (for deuterated drugs)
        # ---------------------------------------------------------------------
        def clean_isotopes_and_hs(m: Chem.Mol) -> Chem.Mol:
            m_copy = Chem.Mol(m)
            for a in m_copy.GetAtoms():
                if a.GetIsotope() != 0:
                    a.SetIsotope(0)
            return Chem.RemoveHs(m_copy)

        notaut_mol = clean_isotopes_and_hs(chosen_frag)
        smi_notaut = Chem.MolToSmiles(notaut_mol, isomericSmiles=True)
        res["parent_smiles_notaut"] = smi_notaut

        # ---------------------------------------------------------------------
        # Step 5: Tautomer canonicalization & comparison
        # ---------------------------------------------------------------------
        taut_mol = te.Canonicalize(notaut_mol)
        taut_mol = clean_isotopes_and_hs(taut_mol)
        smi_taut = Chem.MolToSmiles(taut_mol, isomericSmiles=True)
        res["parent_smiles_taut"] = smi_taut
        res["parent_smiles"] = smi_taut

        if smi_notaut != smi_taut:
            res["tautomer_changed"] = True
            imidic_smarts = Chem.MolFromSmarts("[CX3](-[OX2H])=[NX2]")
            if taut_mol.HasSubstructMatch(imidic_smarts) and not notaut_mol.HasSubstructMatch(imidic_smarts):
                res["imidic_acid_formed"] = True

        p_inchikey = Chem.MolToInchiKey(taut_mol)
        res["parent_inchikey"] = p_inchikey
        res["connectivity_key"] = p_inchikey[:14] if p_inchikey else ""
        res["parent_formula"] = rdMolDescriptors.CalcMolFormula(taut_mol)
        res["parent_mw"] = float(Descriptors.MolWt(taut_mol))
        res["heavy_atoms"] = int(taut_mol.GetNumHeavyAtoms())

    except Exception as e:
        res["parse_ok"] = False
        res["error"] = str(e)

    return res


def main():
    print("=" * 80)
    print("STEP 5: STANDARDIZE STRUCTURES & EXTRACT PARENTS (REFINED)")
    print("=" * 80)

    # 1. Load in-scope resolved structures
    assert RESOLVED_PATH.exists(), f"Missing {RESOLVED_PATH}"
    df_res = pd.read_csv(RESOLVED_PATH)
    in_scope_df = df_res[df_res["in_scope"] == True].copy().reset_index(drop=True)
    print(f"Loaded {len(in_scope_df)} in-scope rows from {RESOLVED_PATH}")
    assert len(in_scope_df) == 1167, f"Expected 1,167 in-scope rows, got {len(in_scope_df)}"

    overrides = load_overrides()
    print(f"Loaded {len(overrides)} parent overrides: {list(overrides.keys())}")

    uncharger = rdMolStandardize.Uncharger()
    te = rdMolStandardize.TautomerEnumerator()
    old_lfc = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)

    records = []
    for idx, row in in_scope_df.iterrows():
        smi = str(row["smiles"]) if pd.notna(row["smiles"]) else ""
        dname = str(row["drug_name"])
        std_info = standardize_molecule(smi, dname, overrides, uncharger, te, old_lfc)

        rec = {
            "ltkb_id": row["ltkb_id"],
            "drug_name": row["drug_name"],
            "source": row["source"],
            "label": row["label"],
            "is_new_v2": row["is_new_v2"],
            "diagnostic_agent": row["diagnostic_agent"],
            "resolution": row["resolution"],
            "original_mol_weight": row["mol_weight"],
            "parse_ok": std_info["parse_ok"],
            "n_fragments_before": std_info["n_fragments_before"],
            "fragment_flag": std_info["fragment_flag"],
            "frag1_smiles": std_info["frag1_smiles"],
            "frag1_heavy_atoms": std_info["frag1_heavy_atoms"],
            "frag2_smiles": std_info["frag2_smiles"],
            "frag2_heavy_atoms": std_info["frag2_heavy_atoms"],
            "parent_needs_review": std_info["parent_needs_review"],
            "parent_override_applied": std_info["parent_override_applied"],
            "parent_smiles_notaut": std_info["parent_smiles_notaut"],
            "parent_smiles_taut": std_info["parent_smiles_taut"],
            "parent_smiles": std_info["parent_smiles"],
            "parent_inchikey": std_info["parent_inchikey"],
            "connectivity_key": std_info["connectivity_key"],
            "parent_formula": std_info["parent_formula"],
            "parent_mw": std_info["parent_mw"],
            "heavy_atoms": std_info["heavy_atoms"],
            "old_parent_smiles": std_info["old_parent_smiles"],
            "old_connectivity_key": std_info["old_connectivity_key"],
            "old_parent_formula": std_info["old_parent_formula"],
            "tautomer_changed": std_info["tautomer_changed"],
            "imidic_acid_formed": std_info["imidic_acid_formed"],
        }
        records.append(rec)

    df_parents = pd.DataFrame(records)
    assert len(df_parents) == 1167, f"Expected 1,167 output rows, got {len(df_parents)}"
    assert df_parents["ltkb_id"].duplicated().sum() == 0, "Duplicate ltkb_id found!"

    # Save data/processed/parents_candidate.csv
    df_parents.to_csv(OUTPUT_PARENTS_PATH, index=False)
    print(f"Saved: {OUTPUT_PARENTS_PATH.resolve()} ({len(df_parents)} rows)")

    # -------------------------------------------------------------------------
    # Topic 6: Parent changes (old chooser vs new chooser)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 6: PARENT CHANGES (CONNECTIVITY KEY DIFFERS OLD vs NEW)")
    print("=" * 80)
    changed_parents = df_parents[df_parents["connectivity_key"] != df_parents["old_connectivity_key"]].copy()
    print(f"Total rows where parent changed: {len(changed_parents)}")
    for _, r in changed_parents.iterrows():
        lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
        print(f"\n  Drug: {r['drug_name']} (Label: {lbl})")
        print(f"    OLD Parent: {r['old_parent_formula']:<14} | Key: {r['old_connectivity_key']} | SMILES: {r['old_parent_smiles']}")
        print(f"    NEW Parent: {r['parent_formula']:<14} | Key: {r['connectivity_key']} | SMILES: {r['parent_smiles']}")

    cols_ch = [
        "ltkb_id", "drug_name", "label",
        "old_connectivity_key", "old_parent_formula", "old_parent_smiles",
        "connectivity_key", "parent_formula", "parent_smiles"
    ]
    changed_parents[cols_ch].to_csv(PROCESSED_DIR / "review5b_parent_changes.csv", index=False)
    print(f"\nSaved: data/processed/review5b_parent_changes.csv ({len(changed_parents)} rows)")

    # -------------------------------------------------------------------------
    # Topic 7: Tautomer Diagnostic
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 7: TAUTOMER DIAGNOSTIC")
    print("=" * 80)
    taut_changed = df_parents[df_parents["tautomer_changed"] == True].copy()
    print(f"Total parents whose SMILES changed upon tautomer canonicalization: {len(taut_changed)}")

    print("\n25 Examples of Tautomer Changes:")
    print(f"{'drug_name':<26} | {'before_tautomer':<36} | {'after_tautomer':<36}")
    print("-" * 105)
    for _, r in taut_changed.head(25).iterrows():
        print(f"{r['drug_name']:<26} | {r['parent_smiles_notaut'][:35]:<36} | {r['parent_smiles_taut'][:35]:<36}")

    imidic_df = df_parents[df_parents["imidic_acid_formed"] == True].copy()
    print(f"\nRows where C=O became C(O)=N (imidic acid / oxime form) ({len(imidic_df)} rows):")
    for _, r in imidic_df.iterrows():
        print(f"  {r['drug_name']:<26} | Before: {r['parent_smiles_notaut']}")
        print(f"  {' '*26} | After : {r['parent_smiles_taut']}")

    # -------------------------------------------------------------------------
    # Topic 5: fragment_flag list (recomputed after neutralization)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 5: FRAGMENT FLAG LIST (recomputed on neutralized organic fragments)")
    print("=" * 80)
    flagged = df_parents[df_parents["fragment_flag"] == True].copy()
    print(f"Total rows with fragment_flag == True: {len(flagged)}\n")
    for _, r in flagged.iterrows():
        lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
        ratio = r["frag1_heavy_atoms"] / r["frag2_heavy_atoms"] if r["frag2_heavy_atoms"] > 0 else 0
        print(f"  {r['drug_name']:<28} | Label: {lbl:<9} | Ratio: {ratio:.2f}")
        print(f"    Frag 1 ({r['frag1_heavy_atoms']:>2} atoms): {r['frag1_smiles']}")
        print(f"    Frag 2 ({r['frag2_heavy_atoms']:>2} atoms): {r['frag2_smiles']}")
    flagged.to_csv(PROCESSED_DIR / "review5_fragment_flag.csv", index=False)
    print(f"Saved: data/processed/review5_fragment_flag.csv ({len(flagged)} rows)")

    # -------------------------------------------------------------------------
    # Topic 8A: Parent MW histogram
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 8A: PARENT MOLECULAR WEIGHT HISTOGRAM")
    print("=" * 80)
    bins = [0, 200, 500, 1000, 2000, float("inf")]
    bin_labels = ["<200", "200-500", "500-1000", "1000-2000", ">2000"]
    df_parents["mw_bin"] = pd.cut(df_parents["parent_mw"], bins=bins, labels=bin_labels, right=False)
    df_parents["label_str"] = df_parents["label"].map({1.0: "1 (Positive)", 0.0: "0 (Negative)"}).fillna("Ambiguous")
    mw_hist = pd.crosstab(df_parents["mw_bin"], df_parents["label_str"], margins=True)
    print(mw_hist.to_string())

    # Heavy parents > 1000
    heavy_parents = df_parents[df_parents["parent_mw"] > 1000].sort_values("parent_mw", ascending=False).copy()
    print(f"\nRows with parent_mw > 1000 ({len(heavy_parents)} rows):")
    print(f"{'drug_name':<30} {'label':>10} {'parent_mw':>12} {'original_mw':>12}")
    print("-" * 68)
    for _, r in heavy_parents.iterrows():
        lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
        print(f"{r['drug_name']:<30} {lbl:>10} {r['parent_mw']:>12.2f} {r['original_mol_weight']:>12.2f}")
    heavy_parents.to_csv(PROCESSED_DIR / "review5_heavy_parents.csv", index=False)
    print(f"Saved: data/processed/review5_heavy_parents.csv ({len(heavy_parents)} rows)")

    # Rows 900 to 1000
    mw_900_1000 = df_parents[(df_parents["parent_mw"] >= 900) & (df_parents["parent_mw"] <= 1000)].sort_values("parent_mw", ascending=False).copy()
    print(f"\nRows with parent_mw between 900 and 1000 ({len(mw_900_1000)} rows):")
    for _, r in mw_900_1000.iterrows():
        lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
        print(f"  {r['drug_name']:<30} | label: {lbl:<9} | parent_mw: {r['parent_mw']:.2f}")
    mw_900_1000.to_csv(PROCESSED_DIR / "review5_parent_mw_900_1000.csv", index=False)
    print(f"Saved: data/processed/review5_parent_mw_900_1000.csv ({len(mw_900_1000)} rows)")

    # -------------------------------------------------------------------------
    # Topic 8B: Duplicate Groups (from new parents)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 8B: DUPLICATE GROUPS (COMPUTED FROM NEW PARENTS)")
    print("=" * 80)
    dup_keys = df_parents[df_parents.duplicated("connectivity_key", keep=False)]["connectivity_key"].unique()
    dup_df = df_parents[df_parents["connectivity_key"].isin(dup_keys)].copy()
    dup_df = dup_df.sort_values(["connectivity_key", "drug_name"])

    groups_conflicting = []
    groups_mixed_ambig = []
    groups_agree = []
    groups_diff_v2 = []

    for k in dup_keys:
        grp = df_parents[df_parents["connectivity_key"] == k]
        labels = grp["label"].dropna().unique()
        has_ambig = grp["label"].isna().any()
        has_labeled = len(labels) > 0
        diff_v2 = len(grp["is_new_v2"].unique()) > 1

        if diff_v2:
            groups_diff_v2.append(k)

        if len(labels) > 1:
            groups_conflicting.append(k)
        elif has_labeled and has_ambig:
            groups_mixed_ambig.append(k)
        elif has_labeled and not has_ambig:
            groups_agree.append(k)
        elif not has_labeled:
            groups_agree.append(k)

    print(f"Total duplicate connectivity groups : {len(dup_keys)}")
    print(f"Total rows involved in duplicate grps: {len(dup_df)}")

    print(f"\n1. Groups with conflicting 0/1 labels ({len(groups_conflicting)} groups):")
    for k in groups_conflicting:
        grp = df_parents[df_parents["connectivity_key"] == k]
        print(f"  Connectivity Key: {k} (Parent: {grp.iloc[0]['parent_formula']}, MW: {grp.iloc[0]['parent_mw']:.1f})")
        for _, r in grp.iterrows():
            lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
            print(f"    - {r['drug_name']:<30} | label: {lbl:<9} | is_new_v2: {str(r['is_new_v2']):<5} | source: {r['source']}")

    print(f"\n2. Groups that mix a labeled and an ambiguous member ({len(groups_mixed_ambig)} groups):")
    for k in groups_mixed_ambig:
        grp = df_parents[df_parents["connectivity_key"] == k]
        print(f"  Connectivity Key: {k}")
        for _, r in grp.iterrows():
            lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
            print(f"    - {r['drug_name']:<30} | label: {lbl:<9} | is_new_v2: {str(r['is_new_v2']):<5} | source: {r['source']}")

    print(f"\n3. Groups where all labeled members agree ({len(groups_agree)} groups):")
    for k in groups_agree:
        grp = df_parents[df_parents["connectivity_key"] == k]
        lbl_rep = str(grp.iloc[0]["label"]).replace(".0", "") if pd.notna(grp.iloc[0]["label"]) else "ambiguous"
        print(f"  Connectivity Key: {k} (Agreed label: {lbl_rep})")
        for _, r in grp.iterrows():
            lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
            print(f"    - {r['drug_name']:<30} | label: {lbl:<9} | is_new_v2: {str(r['is_new_v2']):<5} | source: {r['source']}")

    print(f"\n4. Groups where is_new_v2 differs between members ({len(groups_diff_v2)} groups):")
    for k in groups_diff_v2:
        grp = df_parents[df_parents["connectivity_key"] == k]
        print(f"  Connectivity Key: {k}")
        for _, r in grp.iterrows():
            lbl = str(r["label"]).replace(".0", "") if pd.notna(r["label"]) else "ambiguous"
            print(f"    - {r['drug_name']:<30} | label: {lbl:<9} | is_new_v2: {str(r['is_new_v2']):<5} | source: {r['source']}")

    dup_df.to_csv(PROCESSED_DIR / "review5_duplicate_groups.csv", index=False)
    print(f"\nSaved: data/processed/review5_duplicate_groups.csv ({len(dup_df)} rows)")

    # -------------------------------------------------------------------------
    # Topic 8C: Spot Checks
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 8C: SPOT CHECKS")
    print("=" * 80)
    spot_names = [
        "Abacavir sulfate", "Atorvastatin calcium", "Methotrexate sodium",
        "Valproic acid", "Divalproex sodium", "Pyrvinium pamoate",
        "Deutetrabenazine", "Tetrabenazine", "Omeprazole", "Esomeprazole sodium",
        "Isoniazid", "Rifampin", "Pyrazinamide", "Carbamazepine",
        "Phenytoin", "Ethambutol hydrochloride"
    ]
    spot_rows = []
    print(f"{'drug_name':<26} | {'parent_mw':>9} | {'parent_smiles':<45}")
    print("-" * 88)
    for name in spot_names:
        m = df_parents[df_parents["drug_name"].str.lower() == name.lower()]
        if m.empty:
            print(f"{name:<26} | {'not found':>9} | not found")
            spot_rows.append({
                "drug_name": name, "original_smiles": "not found",
                "parent_smiles": "not found", "parent_mw": None
            })
        else:
            r = m.iloc[0]
            orig_smi = in_scope_df[in_scope_df["ltkb_id"] == r["ltkb_id"]].iloc[0]["smiles"]
            p_mw_str = f"{r['parent_mw']:.2f}" if pd.notna(r['parent_mw']) else "None"
            print(f"{r['drug_name']:<26} | {p_mw_str:>9} | {r['parent_smiles'][:45]}")
            spot_rows.append({
                "drug_name": r["drug_name"],
                "original_smiles": orig_smi,
                "parent_smiles": r["parent_smiles"],
                "parent_mw": r["parent_mw"]
            })
    pd.DataFrame(spot_rows).to_csv(PROCESSED_DIR / "review5_spot_checks.csv", index=False)
    print(f"Saved: data/processed/review5_spot_checks.csv")

    # -------------------------------------------------------------------------
    # Topic 8D: Dry-Run Funnel
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("TOPIC 8D: DRY-RUN FUNNEL (NOTHING APPLIED YET)")
    print("=" * 80)
    total_in = len(df_parents)
    n_parse_fail = len(df_parents[~df_parents["parse_ok"]])

    valid_df = df_parents[df_parents["parse_ok"]].copy()
    heavy_mask = valid_df["parent_mw"] > 1000
    n_heavy = heavy_mask.sum()
    under_1000 = valid_df[~heavy_mask].copy()

    dup_counts = under_1000.groupby("connectivity_key")
    conflict_rows = 0
    collapsed_dup_removed = 0
    retained_rows = []

    for ckey, grp in dup_counts:
        lbls = grp["label"].dropna().unique()
        if len(lbls) > 1:
            conflict_rows += len(grp)
        else:
            collapsed_dup_removed += (len(grp) - 1)
            labeled_members = grp[grp["label"].notna()]
            rep = labeled_members.iloc[0] if not labeled_members.empty else grp.iloc[0]
            retained_rows.append(rep)

    df_retained = pd.DataFrame(retained_rows)

    labeled_retained = df_retained[df_retained["label"].notna()]
    pos_retained = (labeled_retained["label"] == 1).sum()
    neg_retained = (labeled_retained["label"] == 0).sum()
    ambig_retained = df_retained["label"].isna().sum()
    pos_rate = (pos_retained / len(labeled_retained)) * 100 if len(labeled_retained) > 0 else 0.0

    print(f"  Rows in (in-scope candidate pool)       : {total_in}")
    print(f"  Rows with parse failure                 : {n_parse_fail}")
    print(f"  Rows with parent_mw > 1000              : {n_heavy}")
    print(f"  Rows in label-conflict duplicate groups : {conflict_rows}")
    print(f"  Rows removed by collapsing same-label   : {collapsed_dup_removed}")
    print(f"  Remaining rows retained                 : {len(df_retained)}")
    print(f"    - Positives (1)                       : {pos_retained}")
    print(f"    - Negatives (0)                       : {neg_retained}")
    print(f"    - Ambiguous                           : {ambig_retained}")
    print(f"    - Positive rate among labeled drugs   : {pos_rate:.2f}%")

    funnel_data = [
        {"stage": "rows_in", "count": total_in},
        {"stage": "parse_failures", "count": n_parse_fail},
        {"stage": "parent_mw_above_1000", "count": int(n_heavy)},
        {"stage": "label_conflict_rows", "count": conflict_rows},
        {"stage": "collapsible_duplicate_rows_removed", "count": collapsed_dup_removed},
        {"stage": "retained_total", "count": len(df_retained)},
        {"stage": "retained_positives", "count": int(pos_retained)},
        {"stage": "retained_negatives", "count": int(neg_retained)},
        {"stage": "retained_ambiguous", "count": int(ambig_retained)},
        {"stage": "retained_positive_rate_pct", "count": round(pos_rate, 2)},
    ]
    pd.DataFrame(funnel_data).to_csv(PROCESSED_DIR / "review5_dry_run_funnel.csv", index=False)
    print(f"Saved: data/processed/review5_dry_run_funnel.csv")

    print("\n" + "=" * 80)
    print("STEP 5 REFINED STANDARDIZATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
