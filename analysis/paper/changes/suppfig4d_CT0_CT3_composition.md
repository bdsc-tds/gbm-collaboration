# suppfig4d_CT0_CT3_composition.ipynb (new, reconstruction)

No original code found. Implements Methods "Neighborhood cell type composition and marker scoring" (CT_0/CT_3 log2FC vs mean of the other eight neighborhoods, union of top-5 enriched/depleted with ≥ 0.3 % composition; compartment pp CT_0 vs all other neighborhoods).
Data: `AllDataCenters.load_cell2location(sections, abundance="q05_cell_abundance_w_sf")` (identical to the deconv_cohort_level_4 obs of load_spatial_adata), frozen k = 10 CT labels.
Verification: lymphoid + dendritic +7.09 pp (text +7.1). Mesenchymal −6.7 pp not reproduced (MES-like −8.3, MES-like + MES2-like −11.6, malignant −16.2).
