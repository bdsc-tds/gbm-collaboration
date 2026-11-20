#!/usr/bin/env python3
"""
Assemble Supplementary Figure 4 (cell-type and ECM neighborhoods) on two A4 pages (not in the manuscript PDF; letters
from the panel ids; g = the concordance heatmap and its per-sample distribution).
Panels: fig4ab_suppfig4abg_Cell_Neighborhoods, fig4h_suppfig4cd_neighborhood_survival, suppfig4d_CT0_CT3_composition,
fig4e_suppfig4e_colocalization_plots, fig4cdg_suppfig4f_ECM_niches_v4 (symlink_config.yml suppfigure4).
Run (from manuscript/): ./symlink_figures.sh suppfigure4, then from the repo root:
  python manuscript/compose/suppfigure4.py
Output: manuscript/composed/suppfigure4.pdf, suppfigure4_p1..2.png; compose/order/suppfigure4.order.yml.
"""
from _common import cli, compose, flow

ROWS = [
    [("a", "SuppFig4a.pdf")],
    [("b", "SuppFig4b.pdf", 3.0)],
    [("c", "SuppFig4c.pdf"), ("d", "SuppFig4d.pdf")],
    None,
    [("e", "SuppFig4e.pdf")],
    [("f", "SuppFig4f.pdf", 2.5)],
    [("g", "SuppFig4g.pdf", 3.6)],
    [("g", "SuppFig4g_per_sample.pdf", 2.0)],
]

if __name__ == "__main__":
    args = cli("suppfigure4")
    compose("suppfigure4", flow(ROWS, args.panel_dir), args)
