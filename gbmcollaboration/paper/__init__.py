"""Code shared by several paper analyses in analysis/paper/ (Watson et al., MOSAIC GBM).

Data loading uses the regular loaders (gbmcollaboration.AllDataCenters / DataCenter). This
subpackage only holds what more than one analysis needs:

config         paths of frozen inputs, loader flags, constants (exclusion lists live in definitions.py)
ids            patient / block / Visium-section identifiers
cohort         the paper cohort (273 patients, 260 in survival analyses) and covariate coding
annotation     cell-type annotation key
visium         the paper's Visium sections; cell2location merge into obs
til            TIL zones and TIL-high / TIL-low sections
neighborhoods  spot tables and per-neighbourhood sums for the ECM / cell-type neighbourhood analyses
snrna          curated cell-type levels from the annotation key
output         output folders, figure and source-data saving (panels/), package versions
"""
