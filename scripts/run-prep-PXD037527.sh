#!/usr/bin/env bash
data_directory=$(realpath "data")
results_directory=$(realpath "results")

scripts/prepare.sh \
  -r "$results_directory/PXD037527/input-files" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_1.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_2.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_3.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_4.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_5.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_6.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_7.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_8.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_2215_HeYe_1.raw" \
  "$data_directory/PXD037527/Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_2215_HeYe_2.raw"

fast_files=(
  "$data_directory/PXD057758/ConcatenatedHumanEntrapmentProteins_50percent.fasta"
  "$data_directory/PXD057758/uniprot_SCerevisiae_6k_03_2024.fasta"
)

for fasta in "${fast_files[@]}"; do
  base=$(basename "$fasta")

  if [[ ! -e "$results_directory/PXD037527/input-files/$base" ]]; then
    cp "$fasta" "$results_directory/PXD037527/input-files"
  fi
done
