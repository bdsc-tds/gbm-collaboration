#!/usr/bin/env python3
"""
Assemble Supplementary Figure 5 (TIL groups and ECM neighborhoods by TIL zone) on A4 (not in the manuscript PDF; letters
from the panel ids). Panels: fig5c_suppfig5ab_TIL_supplement, fig5de_suppfig5cde_ECM_supplement (symlink_config.yml suppfigure5).
Run (from manuscript/): ./symlink_figures.sh suppfigure5, then from the repo root:
  python manuscript/compose/suppfigure5.py
Output: manuscript/composed/suppfigure5.{pdf,png}; compose/order/suppfigure5.order.yml.
"""
from _common import cli, compose, flow

ROWS = [
    [("a", "SuppFig5a.pdf", 3.8), ("b", "SuppFig5b.pdf")],
    [("c", "SuppFig5c.pdf")],
    [("d", "SuppFig5d.pdf"), ("e", "SuppFig5e.pdf")],
]

if __name__ == "__main__":
    args = cli("suppfigure5")
    compose("suppfigure5", flow(ROWS, args.panel_dir), args)
