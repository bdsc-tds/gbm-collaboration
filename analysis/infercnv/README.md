<div align="center">
<h1>InferCNV single-cell/WES/Visium analyses</h1>
</div>

# Part I - Single-cell data
Script 'infercnv_analyses.py' is used to run InferCNVpy on single-cell data,
using wrappers found in [infercnv_utils.py](./gbmcollaboration/utils/infercnv_utils.py)

Note: this particular line will generate the bins needed for the WES analysis
```python
infercnv_obj.get_bin_positions(
    save_path=RESULTS_DIR / "genomic_windows.csv"
)
```

# Part II - Comparison with WES data
In order to get the same genomics windows (bins) in the WES anlyses as in the single-cell data analysis you need two things:
- run the script 'generate_cnv_wes_matrix.py', part I, that bin CNV data and compute average log2 copy ratio per bin, using the genomic window previously generated as reference for the genomic bins.
- run infercnv_wes_analyses.py to run inferCNVpy visualisation tools
on the WES CNV matrix previously computed.

The script 'analyses_wes_sc_infercnv.py' allows to derive additionnal exploration of the link between single-cell tumors clones and WES CNV clustering. It allows to derive supplementary figure 3.

# Part III - Visium data
Similary type of analyses has been performed on Visium data in the script 'infercnv_visium_analyses.py'.
