# Decisions Log

*Record every change to labels, splits or metrics here, with the date and the reason.*

- **2026-11-01**: Python 3.12 used instead of 3.11 because 3.11 is not installed on this machine; all required libraries support 3.12.
- **2026-11-01**: Defined temporal holdout using drugs with Comment == "New" (added in DILIrank 2.0; 300 total, 169 labeled: 49 positive, 120 negative). To be evaluated once at the end as a separate experiment by retraining on the 813 older labeled drugs (519 positive, 294 negative).
