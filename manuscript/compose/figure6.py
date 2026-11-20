#!/usr/bin/env python3
"""
Assemble Figure 6 (TIL-zone ligand-receptor axes, Xenium) on A4 in the layout of the published figure (manuscript
PDF page 23): boxes read off that page and scaled to A4 by Published; letters at the published positions.
  a        fig6a_suppfig6_NicheNet_circos (its legends are inside the panel, the published right-hand legend area stays empty)
  b, g     fig6bg_xenium_niches, hBrain sections (the sections of the published panels)
  c-f      microscopy, not produced by code: crops of the published page, with the published (white) letters
Run (from manuscript/): ./symlink_figures.sh figure6, then from the repo root:
  python manuscript/compose/figure6.py
Output: manuscript/composed/figure6.{pdf,png}; compose/order/figure6.order.yml.
"""
from _common import Published, cli, compose

P = Published(23, (0.62, 1.10, 7.75, 10.92))

LAYOUT = [
    ("a", P.box(0.80, 1.15, 4.20, 4.10), "Fig6a.pdf"),
    ("b", P.box(5.00, 3.35, 2.60, 1.75), "Fig6b_hBrain.pdf"),
    ("c*", P.box(0.68, 5.32, 3.48, 1.70), P.crop(0.68, 5.32, 3.48, 1.70)),
    ("d*", P.box(4.21, 5.32, 3.51, 1.70), P.crop(4.21, 5.32, 3.51, 1.70)),
    ("e*", P.box(0.68, 7.08, 3.48, 1.65), P.crop(0.68, 7.08, 3.48, 1.65)),
    ("f*", P.box(4.21, 7.08, 3.51, 1.65), P.crop(4.21, 7.08, 3.51, 1.65)),
    ("", P.box(3.15, 8.74, 2.10, 0.22), P.crop(3.15, 8.74, 2.10, 0.22)),
    ("g", P.box(0.70, 8.97, 6.90, 1.95), "Fig6g_hBrain.pdf"),
]

if __name__ == "__main__":
    compose("figure6", LAYOUT, cli("figure6"))
