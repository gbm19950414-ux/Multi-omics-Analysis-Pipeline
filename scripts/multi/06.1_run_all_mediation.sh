# scripts/multi/run_all_mediation.sh
Rscript scripts/multi/06_mediation_analysis_minimal.R \
  results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv \
  M_score_RNA \
  Y_score_RNA \
  results/multiomics/tables/mediation_rna_only

Rscript scripts/multi/06_mediation_analysis_minimal.R \
  results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv \
  M_score_RNA \
  Y_score_lipid \
  results/multiomics/tables/mediation_rna_to_lipid