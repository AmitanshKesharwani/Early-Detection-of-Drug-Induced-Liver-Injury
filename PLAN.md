# Plan

- **Labels**: DILIrank 2.0, positive = vMost + vLess (568), negative = vNo (414), ambiguous excluded (354). Expected total 982.
- **Success criteria**, fixed in advance:
  - The 95 % confidence interval of scaffold‑split AUROC sits clearly above 0.5.
  - The model beats a dummy classifier.
  - The gap between random‑split and scaffold‑split results is reported.
  - AUPRC is compared against the positive rate of about 58 %.
- **Split types and locked test set rule**: stratified locked test set (15‑20 % of drugs), never used for tuning; also repeated stratified CV, repeated scaffold‑group CV, and Butina cluster split.
- **No hyper‑parameter tuning on the locked test set** ever.

## Temporal holdout (pre-registered)
Drugs with Comment == New (added in DILIrank 2.0; 300 total, 169 labeled: 49 positive, 120 negative) form a temporal holdout. Run once, at the end, as a separate experiment: retrain from scratch on the 813 older labeled drugs, choose hyperparameters only by cross-validation within those drugs, then score the 169 new drugs once. Class balance differs (29% vs 64% positive), so AUROC is the headline metric and AUPRC is read against 0.29. A dip is interpreted with caution: new drugs have shorter exposure histories, so no-concern labels may reflect little evidence.

## Scope (pre-registered)
Scope is therapeutic small organic molecules. Exclusion is decided by what a drug is, never by its DILI label. Excluded categories: proteins, antibodies and biologic-scale peptides; oligonucleotides (antisense, siRNA); polymers, resins, mixtures, complexes, diagnostic kits and nanoparticles; diagnostic radiotracers (primary analysis only; a sensitivity analysis will include them); inorganic drugs with no carbon. Metal-containing organic complexes: decide after reviewing the list. After Step 5 (parent extraction), drugs whose parent molecular weight exceeds 1000 Da are excluded, with non-peptide natural products near the cutoff reviewed separately. The expected totals (982 labeled, 568 positive, 414 negative) will shrink; the new positive rate will be computed from data and recorded here before any modeling. All seven target drugs (isoniazid, rifampin, pyrazinamide, methotrexate, valproic acid, carbamazepine, phenytoin) must remain in the retained set; if any is excluded, stop and tell me.

