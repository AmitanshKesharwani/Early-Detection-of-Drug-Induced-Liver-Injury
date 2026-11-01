"""
test_structures.py
==================
Tests for data/processed/structures_raw.csv.

Asserts:
1. Total row count is exactly 1,336 (982 labeled, 354 ambiguous).
2. No duplicate ltkb_id.
3. Preserved label distribution: 568 positives (1) and 414 negatives (0).
4. Ambiguous drugs have null/blank labels.
5. All match_method values are valid.
6. pubchem_failures.csv exactly matches the not_found rows.
"""

from pathlib import Path
import pandas as pd
import pytest

STRUCTURES_PATH = Path("data/processed/structures_raw.csv")
FAILURES_PATH = Path("data/processed/pubchem_failures.csv")


@pytest.fixture(scope="module")
def structures_df():
    assert STRUCTURES_PATH.exists(), f"Missing {STRUCTURES_PATH}"
    df = pd.read_csv(STRUCTURES_PATH)
    return df


def test_row_counts(structures_df):
    """Output must have exactly 1,336 rows: 982 labeled, 354 ambiguous."""
    assert len(structures_df) == 1336, f"Expected 1,336 rows, got {len(structures_df)}"

    labeled_count = (structures_df["source"] == "labeled").sum()
    ambiguous_count = (structures_df["source"] == "ambiguous").sum()

    assert labeled_count == 982, f"Expected 982 labeled, got {labeled_count}"
    assert ambiguous_count == 354, f"Expected 354 ambiguous, got {ambiguous_count}"


def test_no_duplicate_ltkb_id(structures_df):
    """Every drug must have a unique LTKB ID."""
    duplicates = structures_df[structures_df.duplicated("ltkb_id", keep=False)]
    assert duplicates.empty, f"Found duplicated ltkb_id: {duplicates['ltkb_id'].tolist()}"


def test_preserved_labels(structures_df):
    """568 positives and 414 negatives among labeled drugs; ambiguous have no label."""
    labeled_df = structures_df[structures_df["source"] == "labeled"]
    ambiguous_df = structures_df[structures_df["source"] == "ambiguous"]

    positives = (labeled_df["label"] == 1).sum()
    negatives = (labeled_df["label"] == 0).sum()

    assert positives == 568, f"Expected 568 positives, got {positives}"
    assert negatives == 414, f"Expected 414 negatives, got {negatives}"
    assert positives + negatives == 982, "All labeled drugs must be 0 or 1"

    # Ambiguous rows must have null/blank labels
    assert ambiguous_df["label"].isna().all(), "Ambiguous drugs must have null/blank labels"


def test_valid_match_methods(structures_df):
    """All match methods must belong to the allowed set."""
    valid_methods = {"alias", "exact", "isotope_stripped", "salt_stripped", "not_found"}
    actual_methods = set(structures_df["match_method"].unique())
    invalid = actual_methods - valid_methods
    assert not invalid, f"Unexpected match_method values: {invalid}"


def test_failures_consistency(structures_df):
    """pubchem_failures.csv must match the not_found rows in structures_raw.csv."""
    assert FAILURES_PATH.exists(), f"Missing {FAILURES_PATH}"
    failures_df = pd.read_csv(FAILURES_PATH)

    not_found_df = structures_df[structures_df["match_method"] == "not_found"]
    assert len(failures_df) == len(not_found_df), (
        f"Failures count mismatch: {len(failures_df)} vs {len(not_found_df)}"
    )
    assert set(failures_df["ltkb_id"]) == set(not_found_df["ltkb_id"]), (
        "Failures LTKB IDs do not match not_found LTKB IDs"
    )

