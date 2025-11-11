## Full code for the analysis presented in Honours thesis: "Characterising the Tumour Microenvironment of Responders and Non-Responders to Anti-LAG-3-Based Therapy in Advanced Melanoma." by Nabila Zulkapeli.

1. <B>filter_adata.py</b>: merge AnnData with metadata and filter out cores with low melanoma cells, low total cells, remove high TIL cores, and other QC to aid downstream analysis
2. <b>ICI_intro_stats_trials.ipynb</b>: code to make Fig. 1.3
3. <b>umap_matrixplot.ipynb</b>: code to make Fig. 2.2
4. <b>sccoda_ctp.ipynb</b>: code to compare cell type proportions in responders versus non-responders using <i>scCODA</i> (v0.1.9)
5. <b>pertpy_deg.ipynb</b>: code to run differential gene expression analysis in responders versus non-responders using <i>PyDESeq2</i> (v0.5.3)
6. <b>fgsea_region.R</b>: code to run gene set enrichment analysis in peritumour and high-tumour regions using <i>fgsea</i> (v1.32.4)
7. <b>fgsea_figs.ipynb</b>: code to make Fig. 2.6
8. <b>nhood_niche.ipynb</b>: code to perform neighbourhood enrichment analysis per core and response group, and run cluster assignments (extracting flat clusters from dendrograms in clustermaps) for R and NR niches in Table 3.2
9. <b>liana_ccc.ipynb</b>: code to run cell-cell communication inference using <i>CellPhoneDB</i> (v2) and <i>LIANA</i> (v1.6.1)
<br>
<b>Note</b>: All figures were modified in Adobe Illustrator (2024) after exporting from Python. These parameters were used to ensure .pdf files had editable text and vectors:<br>
1. plt.rcParams['pdf.fonttype'] = 42<br>
2. plt.rcParams['ps.fonttype'] = 42<br>
The plots were edited to create multi-panel figures, change colours (when irrelevant to the data, e.g. bar plots), and improve readability.
