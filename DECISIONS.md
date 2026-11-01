# Decisions Log

*Record every change to labels, splits or metrics here, with the date and the reason.*

- **2026-11-01**: Python 3.12 used instead of 3.11 because 3.11 is not installed on this machine; all required libraries support 3.12.
- **2026-11-01**: Defined temporal holdout using drugs with Comment == "New" (added in DILIrank 2.0; 300 total, 169 labeled: 49 positive, 120 negative). To be evaluated once at the end as a separate experiment by retraining on the 813 older labeled drugs (519 positive, 294 negative).
- **2026-11-01**: Pre-registered scope rule defining domain as therapeutic small organic molecules (MW <= 1000 Da after parent extraction). Excluded based on chemical nature regardless of DILI label: proteins, antibodies, biologic-scale peptides, oligonucleotides, polymers/resins/mixtures/complexes/kits/nanoparticles, diagnostic radiotracers (primary analysis), and inorganic compounds with no carbon. Reason: structure review revealed that PubChem not_found entries (109 drugs) are overwhelmingly biologics, oligonucleotides, and macromolecular mixtures; 84 resolved drugs have MW > 1000 Da (mostly peptides/macromolecules); and 20 radiotracers are diagnostic agents rather than therapeutic small molecules. Seven target drugs must remain retained.
- **2026-11-01**: Resolution and scoping rules established:
  - Multi-match rule: each multi-match drug manually reviewed and assigned an explicit CID override based on exact title match, anhydrous form, or primary active stereoisomer recorded in `aliases.csv`.
  - Ivermectin: represented by its primary active component 22,23-dihydroavermectin B1a (CID 6321424, constituting >= 90% of the natural mixture).
  - Drug-carrier complexes: mapped to the active pharmaceutical moiety (e.g. Tetracycline phosphate complex -> Tetracycline CID 54675776; Phentermine resin complex -> Phentermine CID 4771); clinical DILI label inherited from the listed formulation.
  - Metal complexes and inorganics: excluded as a chemistry class regardless of DILI label (27 drugs sit in that category; 8 further metal-containing drugs were assigned to diagnostic_radiotracer by precedence). Organic drugs where metals serve solely as simple salt counterions (Na, K, Ca, Mg, Li) are retained in scope.
  - Deutetrabenazine: kept in scope (therapeutic deuterated small organic molecule, not a diagnostic radiotracer).
  - 1000 Da parent-weight rule: will be applied mechanically in Step 5 with no exceptions (including natural products like temsirolimus and fidaxomicin).
  - Epoetin alfa and Interferon alfa-2b: excluded as macromolecule_biologic_polymer after safety-net review (PubChem returned non-representative model records for these biologics).
  - Verteporfin: categorized as mixture_no_single_structure (regioisomeric benzoporphyrin mixture without a single molecular structure).
  - Colistimethate sodium and Urofollitropin (macromolecule_biologic_polymer) and Estrogens, conjugated (mixture_no_single_structure): excluded after safety-net review because the 1000 Da rule cannot catch biologics for which PubChem returns a partial record, so weight alone is not a safe screen. Estrogens, conjugated is a multi-component mixture where PubChem returned a single component.
- **2026-11-01**: Scope refinements and Step 5 parent extraction rules:
  - Peptides of 5 or more residues (e.g. Eptifibatide, Pentagastrin) are excluded even below 1000 Da as macromolecule_biologic_polymer, while short peptides (4 residues or fewer) and peptidomimetic protease inhibitors remain in scope.
  - Inorganic salts without C-H or C-C bonds (Sodium bicarbonate, Lanthanum carbonate) are excluded under inorganic_or_metal_complex.
  - Perflutren: excluded under mixture_no_single_structure (microsphere formulation; lipid-encapsulated perflutren gas).
  - Dextromethorphan polistirex: mapped by alias to its active moiety dextromethorphan (CID 5360696, formula C18H25NO).
  - Diagnostic dyes and contrast agents (16 drugs) are retained in the dataset but flagged with diagnostic_agent = True for a sensitivity analysis.
  - Step 5 rules: parent extraction uses the largest organic fragment; any drug whose two largest fragments are within a factor of 1.5 in heavy atoms must be flagged for manual review (e.g. pyrvinium pamoate); duplicate structures are collapsed by connectivity.
