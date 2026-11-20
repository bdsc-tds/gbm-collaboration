#!/usr/bin/env python3
"""
Assemble Figure 1 (cohort, modalities, survival) on A4 in the layout of the published figure (manuscript PDF page 5):
boxes (inches, top-left origin) read off that page and scaled to A4 by Published; letters at the published positions.
  a, c, d, e  not produced by code (schematic, circos UMAP, Visium maps, H&E / annotation): crops of the published page
  b           the five KM panels of fig1b_suppfig1_GBM_MS_Summary as a 3 x 2 grid in the published b box; the sixth cell
              is the published OS ridgeline (crop, not produced by code)
Run (from manuscript/): ./symlink_figures.sh figure1, then from the repo root:
  python manuscript/compose/figure1.py
Output: manuscript/composed/figure1.{pdf,png}; compose/order/figure1.order.yml.
"""
from _common import Published, cli, compose, grid

P = Published(5, (0.62, 1.45, 7.70, 11.10))

LAYOUT = [
    ("a", P.box(0.78, 1.50, 3.32, 2.12), P.crop(0.78, 1.50, 3.32, 2.12)),
    *grid("b", P.box(4.22, 1.50, 3.45, 1.85), [
        "Fig1b_mgmt.pdf", "Fig1b_age.pdf", "Fig1b_sex.pdf", "Fig1b_p53.pdf", "Fig1b_center.pdf",
        P.crop(4.40, 1.52, 3.22, 1.05),
    ], ncol=3),
    ("c", P.box(0.62, 3.82, 4.50, 4.58), P.crop(0.62, 3.82, 4.50, 4.58, letter=[(0.74, 3.92, 0.18, 0.22), (5.03, 3.86, 0.10, 0.25)])),
    ("d", P.box(5.12, 3.82, 2.53, 4.41), P.crop(5.12, 3.82, 2.53, 4.41, letter=(5.05, 3.90, 0.18, 0.22))),
    ("e", P.box(0.78, 8.75, 6.52, 2.30), P.crop(0.78, 8.75, 6.52, 2.30)),
]

if __name__ == "__main__":
    compose("figure1", LAYOUT, cli("figure1"))
