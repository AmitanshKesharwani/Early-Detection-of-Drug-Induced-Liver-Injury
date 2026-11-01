"""Tests for label generation and data cleaning integrity."""

from pathlib import Path
import pandas as pd


def test_labels_file_integrity():
    """Verify labels.csv counts, binary values, lack of duplicates, and temporal holdout split."""
    labels_path = Path("data/processed/labels.csv")
    assert labels_path.exists(), "data/processed/labels.csv does not exist. Run src/01_labels.py first."

    df_labels = pd.read_csv(labels_path)

    # 1. Row counts: 568 positives, 414 negatives, 982 total
    positives = (df_labels["label"] == 1).sum()
    negatives = (df_labels["label"] == 0).sum()
    total = len(df_labels)

    assert total == 982, f"Expected 982 total rows in labels.csv, got {total}"
    assert positives == 568, f"Expected 568 positives in labels.csv, got {positives}"
    assert negatives == 414, f"Expected 414 negatives in labels.csv, got {negatives}"

    # 2. Labels are strictly 0 or 1
    unique_labels = set(df_labels["label"].unique())
    assert unique_labels == {0, 1}, f"Expected labels to only contain {0, 1}, got {unique_labels}"

    # 3. No duplicate drug names (case-insensitive, whitespace-stripped)
    drug_names = df_labels["drug_name"].astype(str).str.strip().str.lower()
    dup_names = df_labels[drug_names.duplicated(keep=False)]
    assert dup_names.empty, f"Found duplicate drug names in labels.csv: {dup_names['drug_name'].tolist()}"

    # 4. is_new_v2 column validation & temporal holdout assertions
    assert "is_new_v2" in df_labels.columns, "is_new_v2 column missing from labels.csv"
    # Ensure is_new_v2 matches comment == 'New'
    assert (df_labels["is_new_v2"] == (df_labels["comment"] == "New")).all()

    labeled_new = df_labels[df_labels["is_new_v2"]]
    labeled_older = df_labels[~df_labels["is_new_v2"]]

    # Assertions for New drugs (169 total: 120 label 0, 49 label 1)
    assert len(labeled_new) == 169, f"Expected 169 labeled New drugs, got {len(labeled_new)}"
    assert (labeled_new["label"] == 0).sum() == 120, (
        f"Expected 120 New negatives (0), got {(labeled_new['label'] == 0).sum()}"
    )
    assert (labeled_new["label"] == 1).sum() == 49, (
        f"Expected 49 New positives (1), got {(labeled_new['label'] == 1).sum()}"
    )

    # Assertions for Older drugs (813 total: 294 label 0, 519 label 1)
    assert len(labeled_older) == 813, f"Expected 813 labeled older drugs, got {len(labeled_older)}"
    assert (labeled_older["label"] == 0).sum() == 294, (
        f"Expected 294 older negatives (0), got {(labeled_older['label'] == 0).sum()}"
    )
    assert (labeled_older["label"] == 1).sum() == 519, (
        f"Expected 519 older positives (1), got {(labeled_older['label'] == 1).sum()}"
    )


def test_ambiguous_drugs_file_integrity():
    """Verify ambiguous_drugs.csv counts, is_new_v2 column, and lack of duplicates."""
    ambiguous_path = Path("data/processed/ambiguous_drugs.csv")
    assert ambiguous_path.exists(), (
        "data/processed/ambiguous_drugs.csv does not exist. Run src/01_labels.py first."
    )

    df_ambiguous = pd.read_csv(ambiguous_path)

    # Ambiguous file has exactly 354 rows
    assert len(df_ambiguous) == 354, f"Expected 354 rows in ambiguous_drugs.csv, got {len(df_ambiguous)}"

    # All labels should be blank/null
    assert df_ambiguous["label"].isna().all(), "Ambiguous drugs should not have assigned labels."

    # is_new_v2 column present and matches comment == 'New'
    assert "is_new_v2" in df_ambiguous.columns, "is_new_v2 column missing from ambiguous_drugs.csv"
    assert (df_ambiguous["is_new_v2"] == (df_ambiguous["comment"] == "New")).all()

    # No duplicate drug names
    drug_names = df_ambiguous["drug_name"].astype(str).str.strip().str.lower()
    dup_names = df_ambiguous[drug_names.duplicated(keep=False)]
    assert dup_names.empty, f"Found duplicate drug names in ambiguous_drugs.csv: {dup_names['drug_name'].tolist()}"


def test_no_overlap_between_labeled_and_ambiguous():
    """Gate check: Ensure no drug appears in both labels.csv and ambiguous_drugs.csv."""
    df_labels = pd.read_csv("data/processed/labels.csv")
    df_ambiguous = pd.read_csv("data/processed/ambiguous_drugs.csv")

    labeled_names = set(df_labels["drug_name"].astype(str).str.strip().str.lower())
    ambiguous_names = set(df_ambiguous["drug_name"].astype(str).str.strip().str.lower())

    overlap = labeled_names.intersection(ambiguous_names)
    assert not overlap, f"Found drugs present in both labeled and ambiguous datasets: {overlap}"
