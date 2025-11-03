## Code for analysis presented in the thesis "Identifying biomarkers of response and resistance to anti-LAG-3-based therapy in advanced melanoma." by Nabila Zulkapeli.

1. <B>filter_adata.py</b>: merge AnnData with metadata and filter out cores with low melanoma cells, low total cells, remove high TIL cores, and other QC to aid downstream analysis
2. <b>nhood_per_core</b>: using <B>Squidpy</b> functions, build spatial neighbours per core instead of across all cores and merge and return results as an AnnData object
