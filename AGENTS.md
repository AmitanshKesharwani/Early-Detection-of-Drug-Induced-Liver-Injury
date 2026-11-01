# DILI drug-level model: project rules

Goal: predict DILI concern for a drug from its chemical structure (Objective 1). The user is a beginner: explain every step in plain language and say what to check afterwards.

## Data and labels
- Main labels: DILIrank 2.0. Positive = vMost or vLess concern. Negative = vNo concern. Ambiguous drugs are excluded from training and testing: save them to data/processed/ambiguous_drugs.csv and log the count. Never put them in the locked test set.
- DILIst: use only after removing every drug that appears in DILIrank 2.0 (match by InChIKey connectivity block). It is a limited external check, never extra training data unless I say so.
- Post-2010 drugs added in DILIrank 2.0: use as a temporal holdout only if the file flags them.
- Never edit data/raw. Record download date, URL and file hash in data/raw/README.md.
- Write outputs to data/processed. Cache every API response in data/cache.

## APIs
- PubChem: at most about 5 requests per second, retry with backoff, cache everything, log failures to a file. Read the current PubChem documentation and print one real response before writing a parser.

## Splits and evaluation
- Split by drug, never by row. Keep a stratified locked test set (15 to 20 percent of drugs), touched once at the end.
- Report repeated stratified CV, repeated scaffold-group CV (several seeds, randomized assignment) and a Butina cluster split.
- Metrics: AUROC, AUPRC, MCC, balanced accuracy, Brier score. Report mean and standard deviation. Random seed 42.
- Scaling, column filtering and any feature selection are fitted inside training folds only, in one scikit-learn Pipeline.
- Never tune on the locked test set. Never change labels, splits or metrics without asking me first.

## Habits
- One step per task. Write scripts in src/, not long notebooks.
- After every step print row counts and class balance.
- Add assertions for gate checks (no duplicates, no drug on both sides of a split, expected counts).
- Do not add dependencies without telling me why.
- Do not commit until I say a step passed.

## Stack and layout
Python 3.12, pandas, numpy, scikit-learn, rdkit, lightgbm, shap, requests, matplotlib, pytest.
Folders: data/raw, data/processed, data/cache, src, tests, notebooks, reports.
