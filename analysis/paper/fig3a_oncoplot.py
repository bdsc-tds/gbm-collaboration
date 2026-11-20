"""Fig 3a: oncoplot of the WES driver alterations, and the alteration frequencies of the text.

Produces:
- Fig 3a: oncoplot of the 12 most frequently altered genes (`oncoplot_top12`), with the
  Site code, sex and MGMT tracks.
- The alteration frequencies of the Results ("EGFR amplification (~52%), CDK6 amplification
  (26%), PTEN LOH (21%), TERT promoter mutations (19%), CDK4 amplification (18%), and NTRK2
  amplification (16%)").

Origin: analysis/mutations/mutations_updated_oncoplot.py. Input: the frozen patient x
alteration matrix results/WES/mutations/patient_gene_matrix.csv (pipelines/mutations.py).

What the published panel is (printed below):
- Samples: the mutation-matrix blocks present in the clinical tables loaded with the loader
  exclusions, without the 13 Charite patients of SURVIVAL_EXCLUSION: 251 blocks = 250
  patients (ER_G_290 has blocks a and b, both kept, as it ran). The x label said
  "Samples (n=261)", hard-coded from the number of clinical rows; the text's "WES data from
  261 patients" inherits that number.
- The percentages of the text are the oncoplot's right-hand column: any alteration in the
  gene over the 251 plotted samples. "EGFR amplification ~52%" is EGFR amplified or mutated
  (EGFR_amp alone: 50.6%); "PTEN LOH 21%" is PTEN_lof (deletion log2 < -0.3 or truncating
  SNV/indel; no LOH call is made).
- driver_frequencies_grouped (not in the paper) divides by the 235 samples with at least
  one alteration, so its bars are higher (EGFR_amp 54.0%).

Changes:
- Plot functions and palettes moved verbatim to gbmcollaboration/paper/wes.py.
- Clinical tables: `cohort.clinical_tables(**config.LOADER_FLAGS)` minus SURVIVAL_EXCLUSION
  (the source's per-centre loop, same 261 rows).
- Oncoplot x label: number of plotted samples instead of the hard-coded 261.
- Figures go to the output folder (`output.setup`), .pdf and .png; Fig 3a also to the panels
  folder.
- Dropped: the co-occurrence heatmap (not in the paper), unused imports (seaborn, survival
  and plotting utilities).
- New cell: alteration frequencies with the three denominators (251 plotted blocks, 235
  altered blocks, 250 patients).
"""

# %% Load packages
from pathlib import Path
from typing import Dict, Optional, Sequence

import pandas as pd

from gbmcollaboration import SURVIVAL_EXCLUSION
from gbmcollaboration.paper import cohort, config, output
from gbmcollaboration.paper.ids import patient_of
from gbmcollaboration.paper.wes import (matrix_to_long_events, parse_alt_column, plot_driver_frequencies,
                                        plot_oncoplot)

# =========================
# COLOR MAP
# =========================
CLINICAL_COLOR_MAP = {
    "Site code": {
        "CH": "#0072B2",   # blue
        "CR": "#E69F00",   # orange
        "ER": "#F0E442",   # yellow
        "NA": "#BDBDBD",   # grey (si besoin)
    },
    "Administrative gender": {
        "Male":   "#0072B2",  # blue
        "Female": "#E69F00",  # orange
        "NA": "#BDBDBD",
    },
    "MGMT promoter methylation": {
        "Methylated":   "#0072B2",  # blue
        "Unmethylated": "#E69F00",  # orange
        "NA": "#BDBDBD",
    },
    "Sample origin": {
        "Primary tumour": "#0072B2",
        "Recurrent tumour": "#E69F00",
        "NA": "#BDBDBD",
    },
}

# =========================
# Config
# =========================
RESULTS_DIR = output.setup("06_oncoplot")

SAVE_PLOTS = False
SAVE_FMT = "pdf"                 # "pdf" or "png"
DPI = 300

TOP_GENES_FREQ = 15              # barplot: top genes by any alteration
TOP_GENES_ONCOPRINT = 12         # oncoprint: number of genes shown

# Percentages quoted in the Results, for the check
TEXT_PERCENT = {"EGFR": 52, "CDK6": 26, "PTEN": 21, "TERTp": 19, "CDK4": 18, "NTRK2": 16}


# =========================
# Main entry (call this)
# =========================
def run_all_visualizations(patient_gene_matrix: pd.DataFrame,
                           clinical_df: Optional[pd.DataFrame] = None,
                           clinical_tracks: Optional[Sequence[str]] = None,
                           clinical_color_map: Optional[Dict[str, Dict[object, str]]] = None,
                           outdir: Path = RESULTS_DIR,
                           save: bool = SAVE_PLOTS,
                           fmt: str = SAVE_FMT,
                           dpi: int = DPI) -> Dict[str, object]:
    """
    Create Nature-style cohort visualizations from a boolean patient x alteration matrix.
    Returns a dict with intermediate tables (useful for debugging/reuse).
    """
    events_long = matrix_to_long_events(patient_gene_matrix)

    # 1) Frequencies
    plot_driver_frequencies(events_long, outdir, save=save, fmt=fmt, dpi=dpi, top_n=TOP_GENES_FREQ)

    # 2) Compact oncoprint
    plot_oncoplot(events_long, patient_gene_matrix, clinical_df,
                 outdir, save=save, fmt=fmt, dpi=dpi, top_n_genes=TOP_GENES_ONCOPRINT,
                 clinical_tracks=clinical_tracks,
                 clinical_color_map=clinical_color_map,
                 panel="Fig 3a")

    return {"events_long": events_long}

# %%
# =========================
# Start plotting
# =========================
patient_gene_matrix = pd.read_csv(config.MUTATION_MATRIX_CSV, index_col=0)

# %%
# -- Clinical data of the three data centers (loader exclusions), without the CR samples
# for which survival data is pending
df_concat = cohort.clinical_tables(**config.LOADER_FLAGS).drop(index=SURVIVAL_EXCLUSION)
print(df_concat['Sample collection chronology'].value_counts())
print(f"Found {len(df_concat)} primary cases:", df_concat["center"].value_counts().to_dict())

# %% Start some processing to ensure the oncoplot constains
# only the samples in clinical_df
clinical_df = df_concat.copy()

# Samples used by the oncoplot
samples_onco = patient_gene_matrix.index.astype(str)

# Samples from the clinical_df
samples_clin = clinical_df.index.astype(str)

# Intersection
common_samples = samples_onco.intersection(samples_clin)

print(f"Samples oncoplot: {len(samples_onco)}")
print(f"Samples clinical: {len(samples_clin)}")
print(f"Samples conserved: {len(common_samples)}")
print(f"Samples out: {len(samples_onco) - len(common_samples)}")
print(f"Patients conserved: {patient_of(common_samples).nunique()}")

patient_gene_matrix_filt = patient_gene_matrix.loc[common_samples]
clinical_df_filt = clinical_df.loc[common_samples]

# Need to remove two redundant columns
patient_gene_matrix_filt.drop(columns=["EGFR", "PTEN"], inplace=True)

# %% Fig 3a
outputs = run_all_visualizations(
    patient_gene_matrix_filt,
    clinical_df_filt,
    clinical_tracks=[
        "Site code",
        "Administrative gender",
        "MGMT promoter methylation",
    ],
    clinical_color_map=CLINICAL_COLOR_MAP,
    outdir=RESULTS_DIR,
    save=True
)

# %% Alteration frequencies of the text, with three denominators
# - oncoplot (right-hand column, the text): the plotted blocks
# - driver_frequencies_grouped: the blocks with at least one alteration
# - one block per patient (ER_G_290a dropped)
events_long = outputs["events_long"]
n_plotted = patient_gene_matrix_filt.shape[0]
n_altered = events_long["sample"].nunique()
per_patient = patient_gene_matrix_filt.drop(index=config.DUPLICATE_BLOCK_DROP, errors="ignore")
print(f"Denominators: {n_plotted} plotted blocks, {n_altered} blocks with an alteration, "
      f"{len(per_patient)} patients")

gene_of = {c: parse_alt_column(c)[0] for c in patient_gene_matrix_filt.columns}


def gene_any(m):
    return m.T.groupby(gene_of).any().T


def freq_table(per_gene: bool) -> pd.DataFrame:
    f = gene_any if per_gene else (lambda m: m.astype(bool))
    return pd.DataFrame({
        "n_samples": f(patient_gene_matrix_filt).sum(),
        f"pct_oncoplot (n={n_plotted})": f(patient_gene_matrix_filt).mean() * 100,
        f"pct_driver_freq_plot (n={n_altered} altered)": f(patient_gene_matrix_filt).sum() / n_altered * 100,
        f"pct_per_patient (n={len(per_patient)})": f(per_patient).mean() * 100,
    }).sort_values("n_samples", ascending=False)


freq = freq_table(per_gene=True)
freq["text"] = pd.Series(TEXT_PERCENT)
print("\nAny alteration per gene (the text's % are the oncoplot column):")
print(freq.head(TOP_GENES_ONCOPRINT).round(1).to_string())

by_class = freq_table(per_gene=False)
print("\nPer alteration class:")
print(by_class.loc[["EGFR_amp", "EGFR_mut", "CDK6_amp", "PTEN_lof", "PTEN_mut", "TERTp_mut", "CDK4_amp",
                    "NTRK2_amp"]].round(1).to_string())
output.save_table(RESULTS_DIR, "alteration_frequencies_by_gene", freq.round(2))
output.save_table(RESULTS_DIR, "alteration_frequencies_by_class", by_class.round(2))

# %%
print(output.versions(["numpy", "pandas", "matplotlib"]).to_string(index=False))
