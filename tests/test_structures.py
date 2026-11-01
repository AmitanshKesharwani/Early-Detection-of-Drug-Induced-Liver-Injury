"""
test_structures.py
==================
Tests for data/processed/structures_raw.csv, structures_resolved.csv,
and scope_exclusions.csv.

Asserts:
1. structures_raw.csv and structures_resolved.csv have exactly 1,336 rows (982 labeled, 354 ambiguous).
2. No duplicate ltkb_id in either table.
3. Preserved label distribution in both tables: 568 positives (1) and 414 negatives (0).
4. Every ltkb_id in scope_exclusions.csv exists in structures_resolved.csv with no duplicates.
5. All primary benchmark target drugs are present and marked in_scope == True.
"""

from pathlib import Path
import pandas as pd
import pytest

RAW_PATH = Path("data/processed/structures_raw.csv")
RESOLVED_PATH = Path("data/processed/structures_resolved.csv")
FAILURES_PATH = Path("data/processed/pubchem_failures.csv")
EXCLUSIONS_PATH = Path("data/processed/scope_exclusions.csv")
ALIASES_PATH = Path("data/processed/aliases.csv")


@pytest.fixture(scope="module")
def raw_df():
    assert RAW_PATH.exists(), f"Missing {RAW_PATH}"
    return pd.read_csv(RAW_PATH)


@pytest.fixture(scope="module")
def resolved_df():
    assert RESOLVED_PATH.exists(), f"Missing {RESOLVED_PATH}"
    return pd.read_csv(RESOLVED_PATH)


@pytest.fixture(scope="module")
def exclusions_df():
    assert EXCLUSIONS_PATH.exists(), f"Missing {EXCLUSIONS_PATH}"
    return pd.read_csv(EXCLUSIONS_PATH)


def test_row_counts(raw_df, resolved_df):
    """Both raw and resolved tables must have exactly 1,336 rows (982 labeled, 354 ambiguous)."""
    for name, df in [("raw", raw_df), ("resolved", resolved_df)]:
        assert len(df) == 1336, f"{name}: expected 1,336 rows, got {len(df)}"
        labeled_count = (df["source"] == "labeled").sum()
        ambiguous_count = (df["source"] == "ambiguous").sum()
        assert labeled_count == 982, f"{name}: expected 982 labeled, got {labeled_count}"
        assert ambiguous_count == 354, f"{name}: expected 354 ambiguous, got {ambiguous_count}"


def test_no_duplicate_ltkb_id(raw_df, resolved_df):
    """Every drug must have a unique LTKB ID."""
    for name, df in [("raw", raw_df), ("resolved", resolved_df)]:
        duplicates = df[df.duplicated("ltkb_id", keep=False)]
        assert duplicates.empty, f"{name}: found duplicated ltkb_id: {duplicates['ltkb_id'].tolist()}"


def test_preserved_labels(raw_df, resolved_df):
    """568 positives and 414 negatives among labeled drugs; ambiguous have no label."""
    for name, df in [("raw", raw_df), ("resolved", resolved_df)]:
        labeled_df = df[df["source"] == "labeled"]
        ambiguous_df = df[df["source"] == "ambiguous"]

        positives = (labeled_df["label"] == 1).sum()
        negatives = (labeled_df["label"] == 0).sum()

        assert positives == 568, f"{name}: expected 568 positives, got {positives}"
        assert negatives == 414, f"{name}: expected 414 negatives, got {negatives}"
        assert positives + negatives == 982, f"{name}: labeled sum must equal 982"
        assert ambiguous_df["label"].isna().all(), f"{name}: ambiguous must have null/blank labels"


def test_scope_exclusions_integrity(resolved_df, exclusions_df):
    """Check scope_exclusions.csv uniqueness, reference integrity, and counts."""
    # No duplicate ltkb_id in exclusions
    dup_exclusions = exclusions_df[exclusions_df.duplicated("ltkb_id", keep=False)]
    assert dup_exclusions.empty, f"Duplicates in scope_exclusions.csv: {dup_exclusions['ltkb_id'].tolist()}"

    # Every ltkb_id in scope_exclusions.csv exists in structures_resolved.csv
    resolved_ids = set(resolved_df["ltkb_id"])
    missing_ids = set(exclusions_df["ltkb_id"]) - resolved_ids
    assert not missing_ids, f"Exclusion IDs not in structures_resolved: {missing_ids}"

    # in_scope == False count matches len(exclusions_df)
    excluded_count_in_resolved = (~resolved_df["in_scope"]).sum()
    assert excluded_count_in_resolved == len(exclusions_df), (
        f"Mismatched exclusion count: {excluded_count_in_resolved} vs {len(exclusions_df)}"
    )


def test_target_drugs_in_scope(resolved_df):
    """Target benchmark hepatotoxicity drugs must be retained in scope."""
    target_drugs = [
        "Isoniazid",
        "Rifampin",
        "Pyrazinamide",
        "Methotrexate sodium",
        "Valproic acid",
        "Carbamazepine",
        "Phenytoin",
        "Ethambutol hydrochloride",
    ]
    for target in target_drugs:
        match = resolved_df[resolved_df["drug_name"] == target]
        assert not match.empty, f"Target drug '{target}' missing from structures_resolved.csv"
        is_in_scope = bool(match.iloc[0]["in_scope"])
        assert is_in_scope, f"Target drug '{target}' must be in_scope == True, but got False!"


@pytest.fixture(scope="module")
def aliases_df():
    aliases_path = Path("data/processed/aliases.csv")
    assert aliases_path.exists(), f"Missing {aliases_path}"
    return pd.read_csv(aliases_path)


def test_aliases_csv(aliases_df):
    """aliases.csv must have exactly 4 columns and 13 data rows; every cid_override is an integer;
    Dolasetron mesylate has cid_override == 3033817."""
    assert list(aliases_df.columns) == ["drug_name", "query_name", "cid_override", "note"], (
        f"Unexpected columns: {list(aliases_df.columns)}"
    )
    assert len(aliases_df) == 13, f"Expected 13 rows, got {len(aliases_df)}"
    # Every cid_override must be an integer (not NaN, not float with decimals)
    for idx, row in aliases_df.iterrows():
        val = row["cid_override"]
        assert pd.notna(val), f"Row {idx}: cid_override is NaN"
        assert int(val) == val, f"Row {idx}: cid_override {val!r} is not an integer"
    # Dolasetron mesylate must use 3033817
    dol = aliases_df[aliases_df["drug_name"] == "Dolasetron mesylate"]
    assert not dol.empty, "Dolasetron mesylate not found in aliases.csv"
    assert int(dol.iloc[0]["cid_override"]) == 3033817, (
        f"Dolasetron mesylate cid_override expected 3033817, got {dol.iloc[0]['cid_override']}"
    )


def test_resolved_row_count_and_labels(resolved_df):
    """structures_resolved.csv must have 1,336 rows, 568 positives, 414 negatives."""
    assert len(resolved_df) == 1336, f"Expected 1,336 rows, got {len(resolved_df)}"
    positives = (resolved_df["label"] == 1).sum()
    negatives = (resolved_df["label"] == 0).sum()
    assert positives == 568, f"Expected 568 positives, got {positives}"
    assert negatives == 414, f"Expected 414 negatives, got {negatives}"


def test_biologics_excluded(resolved_df):
    """Excluded biologics, inorganics, and mixtures must have in_scope == False."""
    excluded_targets = [
        "Epoetin alfa",
        "Interferon alfa-2b",
        "Colistimethate sodium",
        "Urofollitropin",
        "Estrogens, conjugated",
        "Eptifibatide",
        "Pentagastrin",
        "Sodium bicarbonate",
        "Lanthanum carbonate",
        "Perflutren",
        "Tetrabamate",
    ]
    for name in excluded_targets:
        match = resolved_df[resolved_df["drug_name"] == name]
        assert not match.empty, f"'{name}' not found in structures_resolved.csv"
        assert not bool(match.iloc[0]["in_scope"]), f"'{name}' should be in_scope == False"


def test_dextromethorphan_polistirex_resolved(resolved_df):
    """Dextromethorphan polistirex must be in_scope == True with formula C18H25NO."""
    match = resolved_df[resolved_df["drug_name"] == "Dextromethorphan polistirex"]
    assert not match.empty, "Dextromethorphan polistirex missing from structures_resolved.csv"
    row = match.iloc[0]
    assert bool(row["in_scope"]) is True, "Dextromethorphan polistirex should be in_scope"
    assert row["formula"] == "C18H25NO", f"Expected formula C18H25NO, got {row['formula']}"


def test_diagnostic_agent_column(resolved_df):
    """diagnostic_agent column exists and is boolean."""
    assert "diagnostic_agent" in resolved_df.columns, "diagnostic_agent column missing"
    assert resolved_df["diagnostic_agent"].dtype == bool, (
        f"diagnostic_agent dtype is {resolved_df['diagnostic_agent'].dtype}, expected bool"
    )
