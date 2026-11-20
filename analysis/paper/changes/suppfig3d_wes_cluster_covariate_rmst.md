# suppfig3d_wes_cluster_covariate_rmst.ipynb (new, reconstruction)

No original code found. Model = Fig 1b / Supp Fig 1c clinical RMST of fig1b_suppfig1_GBM_MS_Summary (one row per subject, τ = 3 y, each covariate adjusted for the others + centre dummies), on the 264 CNV-clustered subjects, without and with a cluster-2 indicator; plus all 287 subjects (cluster 2 = 0 without WES).
`rmst_adjusted` is a compact version of fig1b's `run_rmst` adjusted branch (arm missing -> dropped; covariates missing -> reference level, as fig1b).
Verification: all-subject model reproduces Fig 1b (+7.1 / −5.4 / +3.1 / +3.7 months). With cluster 2: MGMT +6.7 -> +6.3, age −4.9 -> −4.9, sex +3.5 -> +3.2 (all p < 0.011, unchanged significance); p53 +3.0 (p 0.082) -> +3.7 (p 0.033).
