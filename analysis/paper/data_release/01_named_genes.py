"""Genes named in the manuscript (text + figures), with their source.

Figures are raster images in the PDF, so figure genes are taken from the code/tables
that drew them:
  Fig 2b/2c  marker dot plots      -> Glimmer marker sheet (Supp Table 2a), level-4 markers
  Fig 4d/5e, Supp 5c/5d/5e, 6b     -> brain/GBM matrisome panel (Supp Table 3a)
  Fig 6a, Supp 6a/6b circos        -> NicheNet selected axes (ligand / receptor / target)
  Text                             -> curated list below (acronyms that collide with gene
                                      symbols, e.g. MES, LAM, NES, CAR, were removed)
"""

import pandas as pd

import config as C
from common import out

# Curated from the manuscript body (Results/Discussion/Methods). Greek/alias forms
# mapped to HGNC symbols: TGFβ1/2/3 -> TGFB1/2/3, IL-10 -> IL10, IL-15 -> IL15,
# C1Q -> C1QA/B/C, CD8 -> CD8A/CD8B, H3F3A -> H3-3A (+ legacy symbol kept).
TEXT_GENES = """
ACTA2 C1QA C1QB C1QC CCL2 CCN2 CD4 CD8A CD8B CDK4 CDK6 CLU COL1A1 COL1A2 COL3A1
COL4A1 COL4A2 COL5A1 COL5A2 COL6A1 COL6A2 COL6A3 CTNNB1 EGFR EGR1 FGF2 H3C2 H3C3
H3F3A H3-3A IDH1 IDH2 IGF1 IL10 IL15 IL1B IL2RB IL7R LOX LUM MGMT NTRK2 PDGFA PRELP
PTEN RELA RGS2 SNAI1 SRPX2 STAT3 TERT TGFB1 TGFB2 TGFB3 TREM2 TP53
""".split()


def glimmer_markers() -> set:
    df = pd.read_excel(C.GLIMMER_MARKERS_XLSX, sheet_name=0)
    col = [c for c in df.columns if "Level 4 Markers" in str(c)][0]
    genes = set()
    for v in df[col].dropna().astype(str):
        genes |= {g.strip() for g in v.split(",") if g.strip() and g.strip().lower() != "nan"}
    return genes


def main():
    sources: dict[str, set] = {}

    def add(genes, src):
        for g in genes:
            sources.setdefault(g, set()).add(src)

    add(TEXT_GENES, "manuscript_text")
    add(glimmer_markers(), "fig2_marker_dotplot_supp2a")
    add(pd.read_csv(C.ECM_GENES_CSV)["target"], "matrisome_panel_supp3a")
    for p in C.NICHENET_AXES:
        axes = pd.read_csv(p)
        zone = p.parent.name
        for col in ("ligand", "receptor", "target"):
            add(axes[col].dropna().astype(str), f"nichenet_circos_{zone}_{col}")

    df = (
        pd.DataFrame({"gene": list(sources), "sources": [";".join(sorted(s)) for s in sources.values()]})
        .sort_values("gene")
        .reset_index(drop=True)
    )
    df.to_csv(out("genes", "named_genes.csv"), index=False)
    print(f"{len(df)} named genes written")
    print(df.sources.str.split(";").explode().str.replace(r"_(TILplus|BZ1|BZ2)_", "_", regex=True).value_counts())


if __name__ == "__main__":
    main()
