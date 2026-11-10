# Data

No raw data are redistributed. The pipeline in `pipeline.py` needs **no external data**: genetics,
costs and utilities are in `config/params.yaml`.

To recompute the colocalization step (`circadian_dpyd.coloc`) from primary data, obtain:

| Input | Source | Notes |
|---|---|---|
| Hepatic (liver) DPYD eQTL summary statistics | GTEx v8 (gtexportal.org / dbGaP terms apply) | variant-level `snp, beta, se` |
| Variant-level rhythmicity / circadian association data | a human, variant-level source of your choice | **CircaDB provides gene-level rhythmicity (mouse tissues), not variant-level statistics** |
| Clinical outcomes | TCGA-COAD/READ via GDC (n = 382) | outcomes only |

Put files under `data/raw/` (git-ignored) and call:

```python
import pandas as pd
from circadian_dpyd.coloc import coloc_abf
print(coloc_abf(pd.read_csv("data/raw/eqtl.csv"), pd.read_csv("data/raw/rhythm.csv")))
```
The PP4 = 0.83 in the abstract is stored in `config/params.yaml` as an **abstract-reported input**; it is not recomputed by this repository.
