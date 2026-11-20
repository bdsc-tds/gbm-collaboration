# fig3a_0_mutation_calling.py

Origin: `analysis/mutations/mutations.py` (WES patient x alteration matrix, frozen as
results/WES/mutations/patient_gene_matrix.csv, read by fig3a_oncoplot.py).

Changes:
- Header docstring (output, logic, inputs, runtime); change list removed.
- The gene lists are read next to this file (`GENE_LIST_DIR = Path(__file__).parent`; the source
  read them from the working directory, so it only ran from its own folder).
- Matrix written to `output.setup("pipelines_mutations")`, not over the frozen file.
- Dropped the IPython autoreload cell and the oncoplot-style heatmap (`create_oncoplot`,
  superseded by Fig 3a); unused imports removed.

Package helpers: `output` only.

Verification: re-run (17 min 38 s, S3 VCF reads). Not identical to the frozen matrix:
296 vs 299 samples (ER_G_001ult, ER_G_002ult, ER_G_003ult absent from the current S3 listing)
and MET_amp differs for 36 samples; all other 33 columns identical on the 296 common samples.
The run printed "Could not get coordinates for PIK3CA / TP53" (the gene coordinates come from an
online lookup inside mutation_utils), so the gene-coordinate step is not stable over time. The
code edits (paths, header) do not touch the logic; Fig 3a keeps reading the frozen matrix.
