# Raw Data Documentation

This directory contains the original, immutable raw data files for the DILI project. Per project rules, these files must never be edited directly.

---

## 1. Drug Induced Liver Injury Rank (DILIrank 2.0) Dataset

- **Filename**: `Drug Induced Liver Injury Rank (DILIrank 2.0) Dataset  FDA.xlsx`
- **Source URL**: [https://www.fda.gov/science-research/Liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-dataset](https://www.fda.gov/science-research/Liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-dataset)
- **Download Date**: 2026-11-01
- **File Size**: 110,076 bytes
- **SHA-256 Hash**: `1CA1352FF727AF68E68E250EAE2ED775BCA8492335140AC0AFD2233248694993`
- **Sheet Names**:
  - `version 2` (Main dataset — used for labels)
  - `version 1` (Historical reference — do not use)

### Important Notes:
- The workbook has two sheets. **We use only the sheet named `version 2`**.
- In `version 2`, **row 1 is a title banner** (`Drug Induced Liver Injury Rank (DILIrank) Dataset Ver 2.0 | FDA`) and **row 2 holds the actual headers**. Load with `skiprows=1` (or header at row index 1).
- The sheet `version 1` is historical (DILIrank 1.0) and **must not be used for labels**.
- DILIrank concern labels contain spelling variants (vMOST x2, vNo-DILI-Concern x1); normalized by lowercasing; normalized counts verified as 217/351/414/354.

### Sheet Details:
#### `version 2`
- **Total Rows in Sheet**: 1,338 (Row 1: Title banner, Row 2: Header columns, Rows 3–1338: Data)
- **Data Rows**: 1,336
- **Exact Columns**: `LTKBID`, `CompoundName`, `SeverityClass`, `LabelSection`, `vDILI-Concern`, `Comment`
- **Class Distribution (`vDILI-Concern`)**:
  - `vNo-DILI-concern`: 413
  - `Ambiguous-DILI-concern`: 354
  - `vLess-DILI-concern`: 351
  - `vMost-DILI-concern`: 215
  - `vMOST-DILI-concern`: 2
  - `vNo-DILI-Concern`: 1
  - *Total*: 1,336
- **Comment Breakdown**:
  - `Unchanged`: 987
  - `New`: 300
  - `Revised`: 49

#### `version 1` (Historical — Do Not Use)
- **Total Rows in Sheet**: 1,038 (Row 1: Title banner, Row 2: Header columns, Rows 3–1038: Data)
- **Data Rows**: 1,036
- **Exact Columns**: `LTKBID`, `Compound Name`, `Severity Class`, `Label Section`, `vDILIConcern`, `Version`
- **Class Distribution (`vDILIConcern`)**:
  - `No-DILI-Concern`: 312
  - `Less-DILI-Concern`: 278
  - `Ambiguous DILI-concern`: 254
  - `Most-DILI-Concern`: 192
  - *Total*: 1,036

---

## 2. DILIst Dataset

- **Filename**: `DILIst Supplementary Table.xlsx`
- **Source URL**: [https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-severity-and-toxicity-dilist-dataset](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-severity-and-toxicity-dilist-dataset)
- **Download Date**: 2026-11-01
- **File Size**: 49,205 bytes
- **SHA-256 Hash**: `4331EE9D16AE7641488161E4DC2C603C29E06BAA8667DD16E7F3D635366E7E5E`
- **Sheet Names**:
  - `DILIst`

### Important Notes:
- The column names in the raw file have **trailing spaces** (`"DILIst Classification "`, `"Routs of Administration "`) and **must be stripped** when loaded (e.g. `df.columns = df.columns.str.strip()`).
- Per project rules, DILIst is used only as an external check after removing every drug that appears in DILIrank 2.0 (matched by InChIKey connectivity block), never as extra training data unless explicitly requested.

### Sheet Details:
#### `DILIst`
- **Total Rows in Sheet**: 1,280 (Row 1: Header columns, Rows 2–1280: Data)
- **Data Rows**: 1,279
- **Exact Columns**: `DILIST_ID`, `CompoundName`, `DILIst Classification `, `Routs of Administration `
- **Class Distribution (`DILIst Classification `)**:
  - `1` (Positive / DILI Concern): 768
  - `0` (Negative / No DILI Concern): 511
  - *Total*: 1,279
