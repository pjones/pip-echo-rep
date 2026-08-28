#!/usr/bin/env bash
data_directory=$(realpath "data")
results_directory=$(realpath "results")

mzml_files=()
mapfile -t mzml_files <"$results_directory/PXD037527/input-files/mzml_files.lst"

scripts/proteomicslfq.sh \
  -r "$results_directory/PXD037527/proteomicslfq" \
  -f "$results_directory/PXD037527/input-files/ConcatenatedHumanEntrapmentProteins_50percent.fasta" \
  -f "$results_directory/PXD037527/input-files/uniprot_SCerevisiae_6k_03_2024.fasta" \
  "${mzml_files[@]}"
