#!/usr/bin/env bash
results_directory=$(realpath "results")

bash scripts/run-prep-PXD037527.sh

mzml_files=()
mapfile -t mzml_files <"$results_directory/PXD037527/input-files/mzml_files.lst"

scripts/flashlfq.sh \
  -r "$results_directory/PXD037527/flashlfq" \
  -f "$results_directory/PXD037527/input-files/ConcatenatedHumanEntrapmentProteins_50percent.fasta" \
  -f "$results_directory/PXD037527/input-files/uniprot_SCerevisiae_6k_03_2024.fasta" \
  "${mzml_files[@]}"
