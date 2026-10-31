# Plan

- **Labels**: DILIrank 2.0, positive = vMost + vLess (568), negative = vNo (414), ambiguous excluded (354). Expected total 982.
- **Success criteria**, fixed in advance:
  - The 95 % confidence interval of scaffold‑split AUROC sits clearly above 0.5.
  - The model beats a dummy classifier.
  - The gap between random‑split and scaffold‑split results is reported.
  - AUPRC is compared against the positive rate of about 58 %.
- **Split types and locked test set rule**: stratified locked test set (15‑20 % of drugs), never used for tuning; also repeated stratified CV, repeated scaffold‑group CV, and Butina cluster split.
- **No hyper‑parameter tuning on the locked test set** ever.
