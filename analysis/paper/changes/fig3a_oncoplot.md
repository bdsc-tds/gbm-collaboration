# fig3a_oncoplot.py

Origin: `analysis/mutations/mutations_updated_oncoplot.py` (Fig 3a and the alteration frequencies
of the text).

Changes (on the original file; plot functions in the file as written by the author):
- `_maybe_save` writes .pdf and .png with `output.save` (and to the panels folder when `panel` is
  given); `plot_oncoplot` takes `panel` (`"Fig 3a"` from `run_all_visualizations`).
- Oncoplot x label prints the number of plotted samples (was the hard-coded "Samples (n=261)",
  the number of clinical rows; 251 blocks are plotted).
- `RESULTS_DIR = output.setup("fig3a_oncoplot")` instead of figures/oncoplot.
- Prints the number of patients kept (`ids.patient_of`).
- New: alteration frequencies of the Results with three denominators (`TEXT_PERCENT` in Config;
  `gene_any` / `freq_table` with the other functions at the top), tables
  `alteration_frequencies_by_gene/_by_class.csv`; package versions.
- The functions are no longer in `gbmcollaboration/paper/wes.py` (used only here; the module
  was deleted). The co-occurrence heatmap and the original per-centre clinical loading loop
  (`load_tabular` + SURVIVAL_EXCLUSION on CR) are as in the original.

Package helpers: `ids.patient_of` (block -> patient for the printed patient count; shared ID
parser); `output`.

Verification: re-run (54 s). `oncoplot_top12.png`, `driver_frequencies_grouped.png` and both
frequency tables byte-identical to `figures/paper/06_oncoplot/`; 299 matrix samples, 261 clinical,
251 plotted blocks = 250 patients; EGFR 52.2, CDK6 25.9, PTEN 20.7, TERTp 19.1, CDK4 17.9,
NTRK2 16.3 % (text 52/26/21/19/18/16). `driver_cooccurrence_heatmap` is written again (not a panel).
