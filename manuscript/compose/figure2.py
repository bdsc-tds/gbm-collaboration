#!/usr/bin/env python3
"""
Assemble Figure 2 (snRNA-seq atlas and cell-type modules) on A4 in the layout of the published figure (manuscript PDF
page 7): boxes read off that page and scaled to A4 by Published; letters at the published positions.
  a     fig2a_suppfig2_figure2: the UMAP part of Fig2a.pdf (its cell-type legend left out, as published)
  b, c  fig2bc_suppfig2_marker_dotplots
  d, e  fig2de_celltype_modules
Run (from manuscript/): ./symlink_figures.sh figure2, then from the repo root:
  python manuscript/compose/figure2.py
Output: manuscript/composed/figure2.{pdf,png}; compose/order/figure2.order.yml.
"""
from _common import Published, cli, compose

P = Published(7, (0.65, 1.0, 7.65, 10.95), top=0.15)

LAYOUT = [
    ("a", P.box(0.77, 1.20, 3.75, 3.80), "part:Fig2a.pdf#0.05,0.72,5.25,4.05"),
    ("b", P.box(4.57, 1.20, 3.05, 4.20), "Fig2b.pdf"),
    ("c", P.box(0.77, 5.42, 3.05, 5.50), "Fig2c.pdf"),
    ("d", P.box(4.00, 5.42, 3.65, 3.55), "Fig2d.pdf"),
    ("e", P.box(3.85, 9.05, 3.80, 1.85), "Fig2e.pdf"),
]

if __name__ == "__main__":
    compose("figure2", LAYOUT, cli("figure2"))
