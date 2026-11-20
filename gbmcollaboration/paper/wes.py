"""WES helpers shared by the Fig 3 files (06_oncoplot, 07_infercnv_wes_analyses,
08_analyses_wes_sc_infercnv, 09_correlation_cnv_wes_ecm).

- Oncoplot (Fig 3a): the plot functions of analysis/mutations/mutations_updated_oncoplot.py,
  moved verbatim with their palette constants. Changes: `_maybe_save` writes .pdf and .png
  with `output.save` (and to the panels folder when `panel` is given, new argument of
  `plot_oncoplot`); the oncoplot x label prints the number of plotted samples (it was the
  hard-coded "Samples (n=261)").
- `read_h5ad_obs`: the obs table of an .h5ad without reading X (the inferCNV objects are
  16-25 GB on disk; the scripts used them only for their obs).
- `genes_overlapping_segments`: GENCODE genes overlapping genomic segments; replaces the
  pyranges calls of infercnv_wes_analyses.py (pyranges is not in the environment) with the
  same rule as `PyRanges.join` (0-based half-open intervals, any overlap).
- `permanova`: Bray-Curtis PERMANOVA as scikit-bio's `beta_diversity` + `permanova` (not in
  the environment), used for the PERMANOVA numbers of the Methods.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy.spatial.distance import pdist, squareform

from . import output

# =========================================================================== oncoplot
# Recommended alteration colors (colorblind-safe)
ALTERATION_COLORS = {
    "Amplification": "#D55E00", # red-orange
    "Loss-of-Function": "#0072B2",   # blue
    "Mutation": "#009E73",           # green (teal-ish, softer)
    "Altered": "#999999",            # grey
}


ALTERATION_ORDER = [
    "Mutation",
    "Amplification",
    "Loss-of-Function",
    "Altered",
]

MIN_FREQ_TO_LABEL = 0.0          # keep 0 to label all; or e.g. 2.0 to label only >=2%

# Default clinical tracks to show (will be filtered to existing columns)
DEFAULT_CLINICAL_TRACKS = ["Site code"]

# Colorblind-safe "paper-ready" palette (Okabe-Ito)
# We preferentially use blue/orange/yellow for small categorical tracks.
OKABE_ITO = {
    "black":  "#000000",
    "orange": "#E69F00",
    "sky":    "#56B4E9",
    "green":  "#009E73",
    "yellow": "#F0E442",
    "blue":   "#0072B2",
    "verm":   "#D55E00",
    "purple": "#CC79A7",
    "grey":   "#BDBDBD",
}

def _auto_track_palette(n: int) -> list[str]:
    """Return n distinct, colorblind-safe colors (prefer blue/orange/yellow first)."""
    base = [OKABE_ITO["blue"], OKABE_ITO["orange"], OKABE_ITO["yellow"],
            OKABE_ITO["green"], OKABE_ITO["sky"], OKABE_ITO["verm"], OKABE_ITO["purple"], OKABE_ITO["black"]]
    if n <= len(base):
        return base[:n]
    # if more categories than palette, cycle (still acceptable for many centers)
    out = []
    while len(out) < n:
        out.extend(base)
    return out[:n]

def _align_clinical_df(clinical_df: pd.DataFrame, samples_sorted: Sequence[str]) -> pd.DataFrame:
    """Align clinical_df to the sample order used in the oncoplot, with sanity checks."""
    samples_sorted = [str(s) for s in samples_sorted]
    c = clinical_df.copy()

    # 1) Try index first
    c.index = c.index.astype(str)
    match_idx = sum(s in set(c.index) for s in samples_sorted)

    # 2) If poor match, try common sample-id columns
    if match_idx < 0.5 * len(samples_sorted):
        candidate_cols = [col for col in ["sample", "Sample", "sample_id", "Sample ID", "sampleId", "id", "ID"]
                          if col in c.columns]
        best = (match_idx, None)
        for col in candidate_cols:
            idx = c[col].astype(str)
            m = sum(s in set(idx) for s in samples_sorted)
            if m > best[0]:
                best = (m, col)
        if best[1] is not None and best[0] > match_idx:
            c = c.set_index(best[1])
            c.index = c.index.astype(str)
            match_idx = best[0]

    # Reindex to oncoplot order
    c = c.reindex(samples_sorted)

    # Sanity stats (printed once per plot)
    n_missing = int(c.isna().all(axis=1).sum()) if c.shape[1] else 0
    if match_idx < 0.9 * len(samples_sorted):
        print(f"[oncoplot] Clinical alignment: matched {match_idx}/{len(samples_sorted)} samples "
              f"({match_idx/len(samples_sorted):.1%}). Rows fully missing after reindex: {n_missing}. "
              f"Check that clinical_df index (or a sample-id column) matches patient_gene_matrix.index.")
    return c




# =========================
# Helpers
# =========================
def _maybe_save(fig: plt.Figure, outpath: Path, save: bool, dpi: int = 300, panel: Optional[str] = None) -> None:
    if save:
        output.save(outpath.parent, outpath.stem, fig, panel=panel)  # writes .pdf and .png
    plt.show()
    plt.close(fig)


def parse_alt_column(col: str) -> Tuple[str, str]:
    """
    Parse a column name like 'EGFR_amp' or 'TP53_mut' into (gene, event_class).
    event_class in {"Mutation", "Amplification", "Loss-of-Function"} when possible.
    Fallback: ("<raw>", "Altered").
    """
    c = str(col)

    # common suffix patterns
    # mut: mutation (SNV/indel)
    # amp/gain: amplification
    # lof/del/loss: loss-of-function
    m = re.match(r"^(?P<gene>.+?)_(?P<suffix>mut|amp|gain|lof|del|loss)$", c, flags=re.IGNORECASE)
    if m:
        gene = m.group("gene")
        suf = m.group("suffix").lower()
        if suf == "mut":
            evt = "Mutation"
        elif suf in {"amp", "gain"}:
            evt = "Amplification"
        elif suf in {"lof", "del", "loss"}:
            evt = "Loss-of-Function"
        else:
            evt = "Altered"
        return gene, evt

    # handle promoter-style labels like TERTp_mut already covered by suffix regex
    # fallback: try split last underscore
    if "_" in c:
        gene, suf = c.rsplit("_", 1)
        suf = suf.lower()
        if suf in {"mut", "snv", "indel"}:
            return gene, "Mutation"
        if suf in {"amp", "gain"}:
            return gene, "Amplification"
        if suf in {"lof", "del", "loss"}:
            return gene, "Loss-of-Function"

    return c, "Altered"


def matrix_to_long_events(patient_gene_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Convert boolean matrix (patients x alteration-columns) to long table:
    sample, gene, event
    """
    df = patient_gene_matrix.copy()
    df.index = df.index.astype(str)
    df.index.name = "sample"

    # ensure bool
    df = df.applymap(lambda x: bool(x) if pd.notna(x) else False)

    long = df.reset_index().melt(id_vars="sample", var_name="alteration", value_name="present")
    long = long[long["present"]].copy()
    if long.empty:
        raise ValueError("No TRUE alterations found in patient_gene_matrix.")

    parsed = long["alteration"].apply(parse_alt_column)
    long["gene"] = parsed.apply(lambda x: x[0])
    long["event"] = parsed.apply(lambda x: x[1])
    long = long.drop(columns=["present", "alteration"])
    # unique per sample-gene-event
    long = long.drop_duplicates(["sample", "gene", "event"])
    return long


def pick_top_genes(events_long: pd.DataFrame, top_n: int) -> list[str]:
    """
    Top genes by number of samples with any alteration
    """
    g = (
        events_long.drop_duplicates(["sample", "gene"])
        .groupby("gene")["sample"].nunique()
        .sort_values(ascending=False)
    )
    return g.head(top_n).index.tolist()


def make_sample_order(events_long: pd.DataFrame, samples: Sequence[str], genes_priority: Sequence[str]) -> list[str]:
    """
    Sort samples by presence of alterations in a priority list of genes (any event).
    Samples with alterations in earlier genes come first.
    """
    # binary presence per sample per gene
    pres = (
        events_long.drop_duplicates(["sample", "gene"])
        .assign(v=1)
        .pivot_table(index="sample", columns="gene", values="v", aggfunc="max", fill_value=0)
    )
    pres = pres.reindex(index=pd.Index(samples, name="sample")).fillna(0).astype(int)

    # build lexsort keys: earlier genes higher priority
    cols = [g for g in genes_priority if g in pres.columns]
    if not cols:
        return list(samples)

    key = pres[cols].to_numpy()
    # sort descending: convert to lexsort on reversed keys
    order = np.lexsort(key.T[::-1])
    # lexsort gives ascending; we want altered-first -> reverse within each key by sorting on (1 - key)
    key2 = (1 - key)
    order = np.lexsort(key2.T[::-1])
    return list(pres.index[order].astype(str))


# =========================
# Plot 1: Driver frequency landscape
# =========================
def plot_driver_frequencies(events_long: pd.DataFrame,
                           outdir: Path,
                           save: bool = False,
                           fmt: str = "pdf",
                           dpi: int = 300,
                           top_n: int = 15) -> None:
    """
    Driver frequency landscape (grouped bars).

    For each (gene, event_type) we compute:
        %samples with at least one such event in that gene
    and we display the different event types side-by-side (NOT stacked),
    so values cannot exceed 100%.
    """
    if events_long.empty:
        return

    n_total = events_long["sample"].nunique()

    # de-duplicate to avoid counting the same (sample,gene,event) multiple times
    ev = events_long[["sample", "gene", "event"]].drop_duplicates()

    freq = (
        ev
        .groupby(["gene", "event"])["sample"]
        .nunique()
        .reset_index(name="n_samples")
    )
    freq["pct"] = freq["n_samples"] / n_total * 100

    # pick top genes by overall event presence (sum across event types)
    top_genes = (
        freq.groupby("gene")["n_samples"].sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    freq = freq[freq["gene"].isin(top_genes)]

    # consistent event order (reuse existing categories)
    event_order = ["Mutation", "Loss-of-Function", "Amplification", "Altered"]
    freq["event"] = pd.Categorical(freq["event"], categories=event_order, ordered=True)

    pivot = freq.pivot_table(index="gene", columns="event", values="pct", fill_value=0)
    pivot = pivot.loc[list(top_genes)]  # preserve order

    # --- grouped horizontal bars ---
    genes = list(pivot.index.astype(str))
    events = [e for e in event_order if e in pivot.columns]
    vals = pivot[events].to_numpy()

    n_genes = len(genes)
    n_events = len(events)

    fig_h = max(4.5, 0.32 * n_genes)
    fig, ax = plt.subplots(figsize=(8.5, fig_h), dpi=dpi)

    base_y = np.arange(n_genes)
    bar_h = 0.18 if n_events >= 4 else 0.22
    offsets = (np.arange(n_events) - (n_events - 1) / 2) * bar_h

    for j, ev_name in enumerate(events):
        ax.barh(
            base_y + offsets[j],
            vals[:, j],
            height=bar_h * 0.9,
            label=ev_name,
            color=ALTERATION_COLORS[ev_name])

        # annotate each bar (optional)
        for i in range(n_genes):
            v = vals[i, j]
            if (v > 0) and (v >= MIN_FREQ_TO_LABEL):
                ax.text(v + 0.6, base_y[i] + offsets[j], f"{v:.1f}%", va="center", fontsize=8)

    ax.set_yticks(base_y)
    ax.set_yticklabels(genes)
    ax.invert_yaxis()

    ax.set_xlabel("% of samples")
    ax.set_ylabel("")
    ax.legend(title="", loc="lower right", frameon=False)

    plt.tight_layout()
    _maybe_save(fig, outdir / f"driver_frequencies_grouped.{fmt}", save, dpi=dpi)


def plot_oncoplot(events_long: pd.DataFrame,
                  patient_gene_matrix: pd.DataFrame,
                  clinical_df: Optional[pd.DataFrame],
                  outdir: Path,
                  save: bool = False,
                  fmt: str = "pdf",
                  dpi: int = 300,
                  top_n_genes: int = 12,
                  sort_priority_genes: Optional[Sequence[str]] = None,
                  clinical_tracks: Optional[Sequence[str]] = None,
                  clinical_color_map: Optional[Dict[str, Dict[object, str]]] = None,
                  clinical_na_label: str = "NA",
                  panel: Optional[str] = None) -> None:
    """
    Nature-style oncoplot display (close to common oncoplot figures):
      - stacked burden bar on top (counts per sample by event type)
      - oncoprint heatmap (grey background + colored event rectangles)
      - per-gene frequency percentages on the right
      - optional bottom clinical track for 'Site code' (if present)

    IMPORTANT: This function reuses the alteration event classes produced by parse_alt_column()
    from this script (Mutation / Amplification / Loss-of-Function / Altered).
    """
    # ---- choose genes ----
    top_genes = pick_top_genes(events_long, top_n_genes)

    samples = patient_gene_matrix.index.astype(str).tolist()

    # sample sorting (reuse existing default priority)
    if sort_priority_genes is None:
        sort_priority_genes = ["EGFR", "PDGFRA", "NF1", "PTEN", "CDKN2A", "MDM2", "MDM4", "TP53"]
    samples_sorted = make_sample_order(events_long, samples, sort_priority_genes)

    # ---- event order + colors (keep original color intent) ----
    # "Altered" exists as a fallback; we display it as purple to keep it distinct.
    event_order = ["Mutation", "Amplification", "Loss-of-Function", "Altered"]
    event_colors = ALTERATION_COLORS
    # Keep only events that actually occur, but preserve order
    present_events = [e for e in event_order if e in set(events_long["event"].unique())]

    # ---- subset data ----
    sub = events_long[events_long["gene"].isin(top_genes)].copy()
    sub["sample"] = sub["sample"].astype(str)

    # compute stacked burden counts per sample by event (within shown genes)
    burden_counts = (
        sub.groupby(["sample", "event"])
           .size()
           .unstack("event", fill_value=0)
           .reindex(index=samples_sorted, columns=present_events, fill_value=0)
    )

    # compute matrix of events per (gene,row) x (sample,col)
    # store list of events to allow split-rectangles when multiple events occur
    events_map: Dict[Tuple[str, str], list[str]] = {}
    for r in sub.itertuples(index=False):
        key = (str(r.gene), str(r.sample))
        events_map.setdefault(key, [])
        if r.event not in events_map[key]:
            events_map[key].append(r.event)

    # ---- layout ----
    # rows: top bar, oncoprint, optional clinical tracks
    tracks: list[str] = []
    if clinical_df is not None:
        if clinical_tracks is None:
            clinical_tracks = DEFAULT_CLINICAL_TRACKS
        # keep only columns that exist
        tracks = [t for t in clinical_tracks if t in clinical_df.columns]

    nrows = 2 + len(tracks)
    height_ratios = [1.2, 6.0] + [0.7] * len(tracks)

    fig = plt.figure(figsize=(13.5, 0.32 * len(top_genes) + (2.8 + 0.55 * len(tracks))), dpi=dpi)
    LEGEND_X_FIG = 0.92  # figure coordinate x for per-track legends
    gs = fig.add_gridspec(nrows=nrows, ncols=1, height_ratios=height_ratios, hspace=0.15)

    ax_top = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_tracks = [fig.add_subplot(gs[2 + i, 0], sharex=ax_top) for i in range(len(tracks))]
# ---- top stacked burden bar ----
    x = np.arange(len(samples_sorted))
    bottom = np.zeros(len(samples_sorted), dtype=float)
    for evt in present_events:
        vals = burden_counts[evt].to_numpy()
        ax_top.bar(x, vals, bottom=bottom, width=0.9, label=evt, color=event_colors.get(evt, "black"))
        bottom += vals

    ax_top.set_ylabel("# events", fontsize=10)
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.tick_params(axis="y", labelsize=10)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # ---- oncoprint background ----
    n_genes = len(top_genes)
    n_samples = len(samples_sorted)

    ax.set_xlim(-0.5, n_samples - 0.5 + 10)  # extra room for % labels
    ax.set_ylim(n_genes - 0.5, -0.5)

    # draw grey grid as background cells
    bg = np.ones((n_genes, n_samples))
    ax.imshow(bg, aspect="auto", interpolation="nearest",
              cmap=ListedColormap(["#d9d9d9"]), vmin=0, vmax=1)

    # draw event rectangles on top of grey cells
    gene_to_i = {g: i for i, g in enumerate(top_genes)}
    sample_to_j = {s: j for j, s in enumerate(samples_sorted)}

    for (gene, sample), evts in events_map.items():
        if gene not in gene_to_i or sample not in sample_to_j:
            continue
        i = gene_to_i[gene]
        j = sample_to_j[sample]
        # split cell into vertical slices for multiple event types
        evts_sorted = [e for e in event_order if e in evts] + [e for e in evts if e not in event_order]
        k = len(evts_sorted)
        for idx, evt in enumerate(evts_sorted):
            # left edge of slice within cell [-0.5..+0.5]
            x0 = j - 0.5 + idx * (1.0 / k)
            rect = plt.Rectangle((x0, i - 0.5), 1.0 / k, 1.0,
                                 facecolor=event_colors.get(evt, "black"),
                                 edgecolor="white", linewidth=0.35)
            ax.add_patch(rect)

    # ticks/labels
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(top_genes, fontsize=11)
    ax.set_xticks([])  # too many samples usually
    # x ticks are hidden; place xlabel at figure level to avoid overlap with bottom tracks

    # remove spines for clean look
    for sp in ["top", "right", "bottom"]:
        ax.spines[sp].set_visible(False)

    # ---- right-side gene frequency (%) ----
    # percent of samples with any alteration in this gene (within shown genes)
    gene_any = (
        sub.drop_duplicates(["sample", "gene"])
           .assign(v=1)
           .pivot_table(index="gene", columns="sample", values="v", aggfunc="max", fill_value=0)
           .reindex(index=top_genes, columns=samples_sorted, fill_value=0)
    )
    freq_gene = (gene_any.sum(axis=1) / n_samples) * 100

    x_text = n_samples + 0.8
    for i, f in enumerate(freq_gene.to_numpy()):
        if f >= MIN_FREQ_TO_LABEL:
            ax.text(x_text, i, f"{f:.0f}%", va="center", fontsize=10)

    # ---- legend ----
    from matplotlib.patches import Patch
    legend_elems = [Patch(facecolor=event_colors[e], edgecolor="none", label=e) for e in present_events]
    # Place the alterations legend in figure coordinates, aligned to the top bar axis
    pos_top = ax_top.get_position()
    y_top = pos_top.y1
    fig.legend(handles=legend_elems, title="Alterations",
               loc="upper left",
               bbox_to_anchor=(LEGEND_X_FIG, y_top),
               bbox_transform=fig.transFigure,
               frameon=False, fontsize=10, title_fontsize=10)

        # ---- clinical tracks (categorical) ----
    if clinical_df is not None and tracks and ax_tracks:
        # Align to the oncoplot sample order + sanity check
        c = _align_clinical_df(clinical_df, samples_sorted)

        from matplotlib.patches import Patch

        for t_idx, (track, ax_t) in enumerate(zip(tracks, ax_tracks)):
            s_track = c[track] if track in c.columns else pd.Series(index=c.index, dtype=object)

            # Treat as categorical; make missing explicit to avoid misleading big blocks
            s_track = s_track.astype(object).where(~pd.isna(s_track), other=clinical_na_label)
            s_cat = s_track.astype("category")
            cats = list(s_cat.cat.categories)

            # Build category->color map
            user_map = (clinical_color_map or {}).get(track, {})
            pal = _auto_track_palette(len(cats))
            color_map = {}
            for i, cat in enumerate(cats):
                key = cat
                if key in user_map:
                    color_map[key] = user_map[key]
                else:
                    color_map[key] = pal[i]
            # Force NA label to grey unless user explicitly overrides
            if clinical_na_label in cats and clinical_na_label not in user_map:
                color_map[clinical_na_label] = OKABE_ITO["grey"]

            # Stable order
            cats_order = cats
            s_cat = pd.Categorical(s_track, categories=cats_order, ordered=False)
            codes = pd.Series(s_cat, index=c.index).cat.codes.to_numpy()

            # Render a 1-row heat strip
            cmap = ListedColormap([color_map[cat] for cat in cats_order])
            ax_t.imshow(codes[np.newaxis, :], aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=max(1, len(cats_order)-1))
            ax_t.set_yticks([])
            ax_t.set_xticks([])
            ax_t.set_ylabel(str(track), rotation=0, ha="right", va="center", fontsize=8)
            for sp in ax_t.spines.values():
                sp.set_visible(False)

            # Compact legend (cap at 12 categories) — place in *figure* coords, aligned to the track axis
            patches = [Patch(facecolor=color_map[cat], edgecolor="none", label=str(cat)) for cat in cats_order[:12]]
            if patches:
                pos = ax_t.get_position()  # in figure coordinates
                y = (pos.y0 + pos.y1) / 2.0
                fig.legend(handles=patches,
                           title=str(track),
                           loc="center left",
                           bbox_to_anchor=(LEGEND_X_FIG, y),
                           bbox_transform=fig.transFigure,
                           frameon=False,
                           fontsize=8,
                           title_fontsize=10)

    # Global x label placed below all axes
    fig.supxlabel(f"Samples (n={n_samples})", fontsize=12)

    # leave room on the right for legends and % labels
    plt.tight_layout(rect=[0, 0.03, 0.8, 1])
    _maybe_save(fig, outdir / f"oncoplot_top{top_n_genes}.{fmt}", save, dpi=dpi, panel=panel)


# =========================================================================== inferCNV objects
def read_h5ad_obs(path, columns: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """obs of an .h5ad file (optionally only `columns`), without reading X or layers."""
    with h5py.File(path, "r") as f:
        obs = ad.io.read_elem(f["obs"])
    return obs if columns is None else obs[list(columns)]


# =========================================================================== gene annotation
def read_gtf_genes(gtf_path) -> pd.DataFrame:
    """Gene records of a GTF as pyranges.read_gtf gives them (Start = GTF start - 1):
    Chromosome, Start, End, gene_name, gene_id, Strand."""
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                      names=["Chromosome", "Feature", "Start", "End", "Strand", "attributes"],
                      dtype={"Chromosome": str})
    genes = gtf[gtf["Feature"] == "gene"].copy()
    del gtf
    genes["Start"] = genes["Start"] - 1
    genes["gene_name"] = genes["attributes"].str.extract(r'gene_name "([^"]+)"', expand=False)
    genes["gene_id"] = genes["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
    return genes[["Chromosome", "Start", "End", "gene_name", "gene_id", "Strand"]].drop_duplicates()


def genes_overlapping_segments(segments_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    """Join of segments (Chromosome, Start, End, ...) and genes on any overlap of the 0-based
    half-open intervals, as `PyRanges(segments_df).join(PyRanges(genes_df)).df` (gene
    coordinates suffixed `_b`)."""
    hits = segments_df.merge(genes_df, on="Chromosome", suffixes=("", "_b"))
    return hits[(hits["Start_b"] < hits["End"]) & (hits["End_b"] > hits["Start"])]


# =========================================================================== PERMANOVA
def permanova(X, grouping, permutations: int = 999, seed: int = 0) -> pd.Series:
    """PERMANOVA on the Bray-Curtis dissimilarities of the rows of X, as scikit-bio's
    `permanova(beta_diversity("braycurtis", X, ids), grouping, permutations)`.

    The pseudo-F statistic is deterministic; the p value depends on the permutations
    (seeded here, scikit-bio used the global numpy state).
    """
    d2 = squareform(pdist(np.asarray(X, dtype=float), metric="braycurtis")) ** 2
    grouping = np.asarray(grouping)
    n = len(grouping)
    _, codes = np.unique(grouping, return_inverse=True)
    n_groups = codes.max() + 1
    s_t = d2[np.triu_indices(n, 1)].sum() / n

    def pseudo_f(c):
        sizes = np.bincount(c)
        s_w = sum(np.triu(d2[np.ix_(c == k, c == k)], 1).sum() / sizes[k] for k in range(n_groups))
        return ((s_t - s_w) / (n_groups - 1)) / (s_w / (n - n_groups))

    f_obs = pseudo_f(codes)
    rng = np.random.default_rng(seed)
    f_perm = np.array([pseudo_f(rng.permutation(codes)) for _ in range(permutations)])
    p = ((f_perm >= f_obs).sum() + 1) / (permutations + 1)
    return pd.Series({"method name": "PERMANOVA", "test statistic name": "pseudo-F",
                      "sample size": n, "number of groups": n_groups,
                      "test statistic": f_obs, "p-value": p, "number of permutations": permutations})
