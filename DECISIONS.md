# Decisions Log

*Record every change to labels, splits or metrics here, with the date and the reason.*

- **2026-11-01**: Python 3.12 used instead of 3.11 because 3.11 is not installed on this machine; all required libraries support 3.12.
- **2026-11-01**: Defined temporal holdout using drugs with Comment == "New" (added in DILIrank 2.0; 300 total, 169 labeled: 49 positive, 120 negative). To be evaluated once at the end as a separate experiment by retraining on the 813 older labeled drugs (519 positive, 294 negative).
- **2026-11-01**: Pre-registered scope rule defining domain as therapeutic small organic molecules (MW <= 1000 Da after parent extraction). Excluded based on chemical nature regardless of DILI label: proteins, antibodies, biologic-scale peptides, oligonucleotides, polymers/resins/mixtures/complexes/kits/nanoparticles, diagnostic radiotracers (primary analysis), and inorganic compounds with no carbon. Reason: structure review revealed that PubChem not_found entries (109 drugs) are overwhelmingly biologics, oligonucleotides, and macromolecular mixtures; 84 resolved drugs have MW > 1000 Da (mostly peptides/macromolecules); and 20 radiotracers are diagnostic agents rather than therapeutic small molecules. Seven target drugs must remain retained.
