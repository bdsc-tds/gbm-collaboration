#!/usr/bin/env bash
# Re-run every paper analysis and the data release, in manuscript order.
# Light files run one after the other; heavy files (> 60 GB of memory) run one at a time
# under a lock, in parallel with the light queue. Logs: $LOG_DIR/<file>.log, summary in
# $LOG_DIR/run_all.tsv. Upstream steps (fig*_0_*) are not re-run, except the TIL zones.
#
#   bash analysis/paper/run_all.sh            # everything
#   bash analysis/paper/run_all.sh light      # light queue only
#   bash analysis/paper/run_all.sh heavy      # heavy queue only
set -u
cd "$(dirname "$0")"
BIN="${BIN:-$HOME/.conda/envs/gbm_main/bin}"
LOG_DIR="${LOG_DIR:-/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/paper/_logs}"
LOCK="${LOCK:-/tmp/paper_heavy.lock}"
mkdir -p "$LOG_DIR"

LIGHT=(
  table1a_sample_number_per_modality.ipynb
  table1bc_cohort_description.py
  supptables_1bc_pathology_clinical_variables.ipynb
  fig1b_suppfig1_GBM_MS_Summary.ipynb
  fig2de_celltype_modules.ipynb
  fig3a_oncoplot.py
  fig3b_suppfig3ac_infercnv_wes_analyses.py
  fig3c_suppfig3b_wes_profile_os_rmst.ipynb
  fig3d_correlation_cnv_wes_ecm.py
  suppfig3d_wes_cluster_covariate_rmst.ipynb
  suppfig3ef_analyses_wes_sc_infercnv.py
  fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb
  fig4f_pathologist_annotation_analyses.py
  fig4h_suppfig4cd_neighborhood_survival.ipynb
  suppfig4d_CT0_CT3_composition.ipynb
  fig5_0_TIL_zones.ipynb
  fig5a_TIL_zones_map.ipynb
  fig5b_TIL_survival.ipynb
  fig5b_TIL_spline.ipynb
  fig5c_suppfig5ab_TIL_supplement.ipynb
  fig5fg_TIL_zone_GSEA.ipynb
  supptables_4a_til_zone_celltype_stats.ipynb
  supptables_4a_ecm_til_zone_tables.ipynb
  supptables_4b_gsea_TIL_zones.ipynb
  fig6a_suppfig6_NicheNet_circos.ipynb
)
HEAVY=(
  fig2a_suppfig2_figure2.py
  fig2bc_suppfig2_marker_dotplots.ipynb
  fig4cdg_suppfig4f_ECM_niches_v4.ipynb
  fig4e_suppfig4e_colocalization_plots.ipynb
  fig5de_suppfig5cde_ECM_supplement.ipynb
  fig6bg_xenium_niches.ipynb
  data_release/01_named_genes.py
  data_release/02_clinical.py
  data_release/03_wes.py
  data_release/04_singlecell_pseudobulk.py
  data_release/05_visium_pseudobulk.py
  data_release/06_xenium_tiles.py
  data_release/00_cohort_manifest.py
  data_release/07_release_filter.py
  data_release/98_inventory.py
  data_release/99_sanity_checks.py
)

run() {  # run one file, log it, record status and duration
  local f="$1" name start rc
  name=$(basename "${f%.*}")
  start=$(date +%s)
  echo "$(date '+%F %T') START $f" >> "$LOG_DIR/run_all.tsv.log"
  if [[ "$f" == *.ipynb ]]; then
    "$BIN/jupyter" nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=14400 "$f" \
      > "$LOG_DIR/$name.log" 2>&1
  else
    (cd "$(dirname "$f")" && MPLBACKEND=Agg "$BIN/python" "$(basename "$f")") > "$LOG_DIR/$name.log" 2>&1
  fi
  rc=$?
  printf "%s\t%s\t%d\t%ds\n" "$(date '+%F %T')" "$f" "$rc" "$(( $(date +%s) - start ))" >> "$LOG_DIR/run_all.tsv"
}

light() { for f in "${LIGHT[@]}"; do run "$f"; done; }
heavy() { for f in "${HEAVY[@]}"; do flock "$LOCK" bash -c "$(declare -f run); BIN='$BIN' LOG_DIR='$LOG_DIR' run '$f'"; done; }

case "${1:-all}" in
  light) light ;;
  heavy) heavy ;;
  *) light & heavy & wait ;;
esac
echo "done: $LOG_DIR/run_all.tsv"
