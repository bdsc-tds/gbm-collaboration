"""Shared code of the paper analyses in analysis/ (Watson et al., MOSAIC GBM).

config             paths of frozen inputs, loader flags, constants (exclusion lists live in definitions.py)
ids                patient / block / section identifiers
cohort             clinical tables (as the source cells loaded them, and the 273-patient cohort)
visium             Visium expression (batched, gene subset) and cell2location abundances
annotation         cell-type annotation key
til                TIL zones and TIL-high / TIL-low sections
survival           survRM2 RMST, Kaplan-Meier, Cox PH test
clinical_survival  proportional-hazards diagnostics and tau selection (Supp Fig 1a-b)
neighborhood_rmst  RMST sweep over neighbourhood abundances (Fig 4h, Supp Fig 4c)
neighborhoods      ECM / cell-type neighbourhood plotting helpers
snrna              snRNA-seq loading and plotting helpers
wes                oncoplot, WES / CNV helpers
nichenet           NicheNet circos plots (Fig 6a, Supp Fig 6)
xenium             Xenium tiling and TIL correlation helpers (Fig 6b, 6g)
output             output folders, panel copies, figure saving, package versions
"""
