#!/usr/bin/env bash
data_directory=$(realpath "data")
results_directory=$(realpath "results")

scripts/flashlfq.sh \
  -r "$results_directory/PXD037527/flashlfq" \
  -f "$data_directory/PXD037527/uniprot_reviewed_combined_UP5640_UP2311_2022_08_11.fasta" \
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
