#!/usr/bin/env python3
"""
Assemble Supplementary Figure 2 (snRNA-seq annotation by compartment, levels 1-4) on six A4 pages (not in the manuscript
PDF; letters from the panel ids): 1 tumour UMAPs a and composition bars b beside the marker dot plot c; 2 immune d, e;
3 immune dot plots f; 4 stromal / neuronal g, h; 5-6 their dot plots i (Glia, Neuron, Vascular; Stromal).
Panels: fig2a_suppfig2_figure2, fig2bc_suppfig2_marker_dotplots (symlink_config.yml suppfigure2).
Run (from manuscript/): ./symlink_figures.sh suppfigure2, then from the repo root:
  python manuscript/compose/suppfigure2.py
Output: manuscript/composed/suppfigure2.pdf, suppfigure2_p1..6.png; compose/order/suppfigure2.order.yml.
"""
from _common import cli, compose, flow


def levels(letter, prefix):
    return [(letter, f"{prefix}_level{i}.pdf") for i in (1, 2, 3, 4)]


def umaps_and_bars(umap, bar, pair_bars=True):
    u, b = levels(*umap), levels(*bar)
    bars = [[b[0], b[1]], [(*b[2], 1.5)], [b[3]]] if pair_bars else [[b[0]], [b[1]], [b[2]], [b[3]]]
    return [[(*u[0], 2.3), u[1]], [(*u[2], 2.3), u[3]], *bars]


def pages(panel_dir):
    tumour = flow(umaps_and_bars(("a", "SuppFig2a"), ("b", "SuppFig2b"), pair_bars=False), panel_dir, right=5.15, max_h=2.3)[0]
    tumour += flow([[("c", "SuppFig2c_Malignant.pdf", 11.0)]], panel_dir, left=5.25)[0]
    return [tumour] + flow([
        *umaps_and_bars(("d", "SuppFig2d"), ("e", "SuppFig2e")), None,
        [("f", "SuppFig2f_Lymphoid.pdf", 11.0), ("f", "SuppFig2f_Myeloid.pdf")], None,
        *umaps_and_bars(("g", "SuppFig2g"), ("h", "SuppFig2h")), None,
        [("i", "SuppFig2i_Glia.pdf", 11.0), ("i", "SuppFig2i_Neuron.pdf"), ("i", "SuppFig2i_Vascular.pdf")], None,
        [("i", "SuppFig2i_Stromal.pdf", 11.0)],
    ], panel_dir, max_h=2.6)


if __name__ == "__main__":
    args = cli("suppfigure2")
    compose("suppfigure2", pages(args.panel_dir), args)
