"""Screen drug names for structural/curation flags prior to PubChem querying.

Reads:
- data/processed/labels.csv
- data/processed/ambiguous_drugs.csv

Flags:
- combination: contains '/', ' and ', ' with ', ' + ' or ','
- has_parentheses: contains '(' or ')'
- has_digits: contains any numeric digit
- non_ascii: contains any character with ord > 127
- extra_whitespace: leading, trailing or double spaces
- salt_or_hydrate: contains salt/counterion/hydrate forms
- likely_biologic: names ending with -mab, -cept, -tide, -ase, -kin or similar

Outputs:
- data/processed/name_review.csv
- data/processed/aliases.csv (empty template: drug_name, query_name, note)
"""

from pathlib import Path
import re
import pandas as pd

# Salt / hydrate roots to flag
SALT_HYDRATE_WORDS = [
    "sulfate",
    "hydrochloride",
    "sodium",
    "potassium",
    "calcium",
    "maleate",
    "mesylate",
    "tartrate",
    "acetate",
    "phosphate",
    "citrate",
    "besylate",
    "bromide",
    "chloride",
    "succinate",
    "fumarate",
    "hydrobromide",
    "hydrate",
]

# Regex matches exact words or prefixed variations like disodium, dihydrochloride, monohydrate
SALT_PATTERN = re.compile(
    r"\b\w*(?:" + "|".join(SALT_HYDRATE_WORDS) + r")\b",
    re.IGNORECASE,
)

# Biologic suffix pattern (mab, cept, tide, ase, kin)
BIOLOGIC_PATTERN = re.compile(
    r"\b\w*(?:mab|cept|tide|ase|kin)\b",
    re.IGNORECASE,
)

# Combination delimiters: '/', ' and ', ' with ', ' + ', or ','
COMBINATION_PATTERN = re.compile(
    r"/|\band\b|\bwith\b|\+|\,",
    re.IGNORECASE,
)


def flag_drug_names(df: pd.DataFrame) -> pd.DataFrame:
    """Compute boolean flags for drug names without modifying the names."""
    names = df["drug_name"].astype(str)

    df["combination"] = names.apply(lambda x: bool(COMBINATION_PATTERN.search(x)))
    df["has_parentheses"] = names.apply(lambda x: "(" in x or ")" in x)
    df["has_digits"] = names.apply(lambda x: bool(re.search(r"\d", x)))
    df["non_ascii"] = names.apply(lambda x: any(ord(c) > 127 for c in x))
    df["extra_whitespace"] = names.apply(lambda x: x != x.strip() or "  " in x)
    df["salt_or_hydrate"] = names.apply(lambda x: bool(SALT_PATTERN.search(x)))
    df["likely_biologic"] = names.apply(lambda x: bool(BIOLOGIC_PATTERN.search(x)))

    return df


def main():
    processed_dir = Path("data/processed")
    labels_file = processed_dir / "labels.csv"
    ambiguous_file = processed_dir / "ambiguous_drugs.csv"

    if not labels_file.exists() or not ambiguous_file.exists():
        raise FileNotFoundError("Missing labels.csv or ambiguous_drugs.csv in data/processed.")

    df_labels = pd.read_csv(labels_file)
    df_ambiguous = pd.read_csv(ambiguous_file)

    df_labels["source"] = "labeled"
    df_ambiguous["source"] = "ambiguous"

    df_combined = pd.concat([df_labels, df_ambiguous], ignore_index=True)
    total_drugs = len(df_combined)
    print(f"Loaded {len(df_labels)} labeled drugs and {len(df_ambiguous)} ambiguous drugs.")
    print(f"Total drugs to review: {total_drugs}\n")

    # Add flags
    df_reviewed = flag_drug_names(df_combined)

    flag_columns = [
        "combination",
        "has_parentheses",
        "has_digits",
        "non_ascii",
        "extra_whitespace",
        "salt_or_hydrate",
        "likely_biologic",
    ]

    # Save to data/processed/name_review.csv
    review_output_cols = [
        "ltkb_id",
        "drug_name",
        "source",
        "label",
        "is_new_v2",
        "comment",
    ] + flag_columns

    review_file = processed_dir / "name_review.csv"
    df_reviewed[review_output_cols].to_csv(review_file, index=False)
    print(f"Saved drug name review to: {review_file}\n")

    # Create empty aliases.csv if it doesn't already exist or as specified
    aliases_file = processed_dir / "aliases.csv"
    df_aliases_template = pd.DataFrame(columns=["drug_name", "query_name", "note"])
    df_aliases_template.to_csv(aliases_file, index=False)
    print(f"Created template alias mapping at: {aliases_file}\n")

    # Print summary and examples for each flag
    print("=" * 70)
    print("DRUG NAME REVIEW: SUMMARY & EXAMPLES")
    print("=" * 70)

    for flag in flag_columns:
        matching = df_reviewed[df_reviewed[flag]]
        count = len(matching)
        pct = (count / total_drugs) * 100
        print(f"\n--- Flag: {flag} ---")
        print(f"Count: {count} / {total_drugs} ({pct:.2f}%)")

        if count == 0:
            print("  (None found)")
        else:
            examples = matching["drug_name"].head(15).tolist()
            for i, name in enumerate(examples, 1):
                source = matching.iloc[i - 1]["source"]
                print(f"  {i:2d}. {name} [{source}]")
            if count > 15:
                print(f"  ... and {count - 15} more.")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()

