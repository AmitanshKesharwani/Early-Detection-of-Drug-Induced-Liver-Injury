"""
tests/test_parents.py
=====================
Tests for chemical standardization and data/processed/parents_candidate.csv.
"""

from pathlib import Path
import pandas as pd
import pytest

RESOLVED_PATH = Path("data/processed/structures_resolved.csv")
PARENTS_PATH = Path("data/processed/parents_candidate.csv")
FINAL_PARENTS_PATH = Path("data/processed/parents_final.csv")
HEAVY_REVIEW_PATH = Path("data/processed/review6_heavy_dropped.csv")
COLLAPSED_REVIEW_PATH = Path("data/processed/review6_collapsed.csv")


@pytest.fixture(scope="module")
def full_resolved_df():
    assert RESOLVED_PATH.exists(), f"Missing {RESOLVED_PATH}"
    return pd.read_csv(RESOLVED_PATH)


@pytest.fixture(scope="module")
def in_scope_resolved_df(full_resolved_df):
    return full_resolved_df[full_resolved_df["in_scope"] == True].copy().reset_index(drop=True)


@pytest.fixture(scope="module")
def parents_df():
    assert PARENTS_PATH.exists(), f"Missing {PARENTS_PATH}"
    return pd.read_csv(PARENTS_PATH)


@pytest.fixture(scope="module")
def final_parents_df():
    assert FINAL_PARENTS_PATH.exists(), f"Missing {FINAL_PARENTS_PATH}"
    return pd.read_csv(FINAL_PARENTS_PATH)


def test_standardize_key_targets_mw(parents_df):
    """Abacavir, Atorvastatin, Methotrexate, and Valproic acid parent MW within 0.5."""
    targets = [
        ("Abacavir sulfate", 286.3, 0.5),
        ("Atorvastatin calcium", 558.6, 0.5),
        ("Methotrexate sodium", 454.4, 0.5),
        ("Valproic acid", 144.2, 0.5),
    ]
    for name, expected_mw, tolerance in targets:
        match = parents_df[parents_df["drug_name"] == name]
        assert not match.empty, f"'{name}' missing from parents_candidate.csv"
        p_mw = match.iloc[0]["parent_mw"]
        assert abs(p_mw - expected_mw) <= tolerance, (
            f"{name} parent MW expected {expected_mw} +/- {tolerance}, got {p_mw:.2f}"
        )


def test_tafamidis_iothalamate_formulas(parents_df):
    """Tafamidis meglumine parent formula C14H7Cl2NO3; Iothalamate meglumine parent formula C11H9I3N2O4."""
    taf = parents_df[parents_df["drug_name"] == "Tafamidis meglumine"]
    assert not taf.empty, "Tafamidis meglumine missing"
    assert taf.iloc[0]["parent_formula"] == "C14H7Cl2NO3", (
        f"Expected C14H7Cl2NO3, got {taf.iloc[0]['parent_formula']}"
    )

    iot = parents_df[parents_df["drug_name"] == "Iothalamate meglumine"]
    assert not iot.empty, "Iothalamate meglumine missing"
    assert iot.iloc[0]["parent_formula"] == "C11H9I3N2O4", (
        f"Expected C11H9I3N2O4, got {iot.iloc[0]['parent_formula']}"
    )


def test_override_parent_formulas(parents_df):
    """Fosfomycin tromethamine parent formula C3H7O4P; Methenamine hippurate parent formula C6H12N4."""
    fos = parents_df[parents_df["drug_name"] == "Fosfomycin tromethamine"]
    assert not fos.empty, "Fosfomycin tromethamine missing"
    assert fos.iloc[0]["parent_formula"] == "C3H7O4P", (
        f"Expected C3H7O4P, got {fos.iloc[0]['parent_formula']}"
    )

    meth = parents_df[parents_df["drug_name"] == "Methenamine hippurate"]
    assert not meth.empty, "Methenamine hippurate missing"
    assert meth.iloc[0]["parent_formula"] == "C6H12N4", (
        f"Expected C6H12N4, got {meth.iloc[0]['parent_formula']}"
    )


def test_dimenhydrinate_formula(parents_df):
    """Dimenhydrinate parent formula C17H21NO."""
    dim = parents_df[parents_df["drug_name"] == "Dimenhydrinate"]
    assert not dim.empty, "Dimenhydrinate missing"
    assert dim.iloc[0]["parent_formula"] == "C17H21NO", (
        f"Expected C17H21NO, got {dim.iloc[0]['parent_formula']}"
    )


def test_divalproex_parent_equals_valproic_acid(parents_df):
    """Divalproex sodium parent SMILES equals Valproic acid's."""
    div = parents_df[parents_df["drug_name"] == "Divalproex sodium"]
    val = parents_df[parents_df["drug_name"] == "Valproic acid"]
    assert not div.empty, "Divalproex sodium missing"
    assert not val.empty, "Valproic acid missing"
    assert div.iloc[0]["parent_smiles"] == val.iloc[0]["parent_smiles"], (
        f"Divalproex parent {div.iloc[0]['parent_smiles']} != Valproic acid {val.iloc[0]['parent_smiles']}"
    )


def test_deutetrabenazine_no_explicit_hs_and_matches_tetrabenazine(parents_df):
    """Deutetrabenazine parent contains no explicit [H] atoms and has the same connectivity key as Tetrabenazine."""
    deute = parents_df[parents_df["drug_name"] == "Deutetrabenazine"]
    tetra = parents_df[parents_df["drug_name"] == "Tetrabenazine"]
    assert not deute.empty, "Deutetrabenazine missing"
    assert not tetra.empty, "Tetrabenazine missing"

    deute_smi = deute.iloc[0]["parent_smiles"]
    assert "[H]" not in deute_smi, f"Deutetrabenazine parent still contains [H]: {deute_smi}"

    deute_ckey = deute.iloc[0]["connectivity_key"]
    tetra_ckey = tetra.iloc[0]["connectivity_key"]
    assert deute_ckey == tetra_ckey, (
        f"Connectivity key mismatch: {deute_ckey} vs {tetra_ckey}"
    )


def test_tetrabamate_excluded(full_resolved_df):
    """Tetrabamate is in_scope False."""
    tetra = full_resolved_df[full_resolved_df["drug_name"] == "Tetrabamate"]
    assert not tetra.empty, "Tetrabamate missing from structures_resolved.csv"
    assert not bool(tetra.iloc[0]["in_scope"]), "Tetrabamate must be in_scope == False"


def test_parents_candidate_integrity(in_scope_resolved_df, parents_df):
    """parents_candidate.csv has one row per in-scope drug with no duplicate ltkb_id."""
    assert len(parents_df) == len(in_scope_resolved_df), (
        f"Mismatch: parents_candidate has {len(parents_df)} rows, in_scope has {len(in_scope_resolved_df)}"
    )

    dups = parents_df[parents_df.duplicated("ltkb_id", keep=False)]
    assert dups.empty, f"Found duplicated ltkb_id: {dups['ltkb_id'].tolist()}"

    merged = parents_df.merge(
        in_scope_resolved_df[["ltkb_id", "label", "is_new_v2"]],
        on="ltkb_id",
        suffixes=("_cand", "_res"),
    )
    assert len(merged) == len(in_scope_resolved_df)

    label_match = (
        (merged["label_cand"] == merged["label_res"])
        | (merged["label_cand"].isna() & merged["label_res"].isna())
    )
    assert label_match.all(), "Labels in parents_candidate.csv do not match structures_resolved.csv"

    v2_match = merged["is_new_v2_cand"] == merged["is_new_v2_res"]
    assert v2_match.all(), "is_new_v2 in parents_candidate.csv does not match structures_resolved.csv"


def test_parents_final_counts_and_integrity(final_parents_df, parents_df):
    """Step 6 cutoffs and collapses yield exactly 1,094 rows with expected counts."""
    assert len(final_parents_df) == 1094, f"Expected 1094 rows, got {len(final_parents_df)}"

    # Gate checks
    assert final_parents_df["ltkb_id"].nunique() == len(final_parents_df)
    assert final_parents_df["connectivity_key"].notna().all()

    # Label counts
    pos = (final_parents_df["label"] == 1).sum()
    neg = (final_parents_df["label"] == 0).sum()
    ambig = final_parents_df["label"].isna().sum()
    assert pos == 511, f"Expected 511 positives, got {pos}"
    assert neg == 276, f"Expected 276 negatives, got {neg}"
    assert ambig == 307, f"Expected 307 ambiguous, got {ambig}"


def test_parents_final_mw_cutoff(final_parents_df):
    """No drug in parents_final.csv has parent_mw > 1000."""
    heavy = final_parents_df[final_parents_df["parent_mw"] > 1000]
    assert heavy.empty, f"Found drugs with parent_mw > 1000: {heavy['drug_name'].tolist()}"


def test_conflict_group_both_retained(final_parents_df):
    """Amphetamine sulfate (pos) and Dextroamphetamine sulfate (neg) are both retained."""
    amp = final_parents_df[final_parents_df["drug_name"] == "Amphetamine sulfate"]
    damp = final_parents_df[final_parents_df["drug_name"] == "Dextroamphetamine sulfate"]
    assert not amp.empty, "Amphetamine sulfate missing from parents_final.csv"
    assert not damp.empty, "Dextroamphetamine sulfate missing from parents_final.csv"
    assert amp.iloc[0]["label"] == 1.0
    assert damp.iloc[0]["label"] == 0.0
    assert amp.iloc[0]["connectivity_key"] == damp.iloc[0]["connectivity_key"]


def test_same_label_collapse_example(final_parents_df):
    """Valproic acid is kept, Divalproex sodium is collapsed, and collapsed_from records both."""
    val = final_parents_df[final_parents_df["drug_name"] == "Valproic acid"]
    div = final_parents_df[final_parents_df["drug_name"] == "Divalproex sodium"]
    assert not val.empty, "Valproic acid missing from parents_final.csv"
    assert div.empty, "Divalproex sodium should have been collapsed into Valproic acid"
    assert "Divalproex sodium" in str(val.iloc[0]["collapsed_from"])

