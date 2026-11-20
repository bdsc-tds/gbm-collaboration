#!/usr/bin/env python3
"""
Assemble Figure 5 (TIL zones) on A4 in the layout of the published figure (manuscript PDF page 20): boxes read off
that page and scaled to A4 by Published; letters at the published positions.
  a  fig5a_TIL_zones_map
  b  fig5b_TIL_survival
  c  fig5c_suppfig5ab_TIL_supplement
  d, e  fig5de_suppfig5cde_ECM_supplement
  f, g  fig5fg_TIL_zone_GSEA (g: TIL+ vs BZ1 above BZ1 vs BZ2)
Run (from manuscript/): ./symlink_figures.sh figure5, then from the repo root:
  python manuscript/compose/figure5.py
Output: manuscript/composed/figure5.{pdf,png}; compose/order/figure5.order.yml.
"""
from _common import LETTER_W, Published, cli, compose

P = Published(20, (0.55, 1.05, 7.60, 11.0))

LAYOUT = [
    ("a", P.box(0.62, 1.10, 3.45, 2.15), "Fig5a.pdf"),
    ("b", P.box(4.22, 1.10, 3.30, 2.15), "Fig5b.pdf"),
    ("c", P.box(0.62, 3.38, 6.95, 1.80), "Fig5c.pdf"),
    ("d", P.box(0.62, 5.25, 4.95, 1.75), "Fig5d.pdf"),
    ("e", P.box(5.62, 5.25, 1.90, 1.75), "Fig5e.pdf"),
    ("f", P.box(0.60, 7.20, 3.30, 3.80), "Fig5f.pdf"),
    ("g", P.box(3.93, 7.20, 3.62, 1.85), "Fig5g_TIL+_vs_BZ1.pdf"),
    ("g", P.box(3.93 + LETTER_W / P.k, 9.10, 3.62 - LETTER_W / P.k, 1.85), "Fig5g_BZ1_vs_BZ2.pdf"),
]

if __name__ == "__main__":
    compose("figure5", LAYOUT, cli("figure5"))
