#!/usr/bin/env bash

data_directory=$(realpath "data")
results_directory=$(realpath "results")

mzml_files=()
mapfile -t mzml_files <"$results_directory/PXD037527/input-files/mzml_files.lst"

bash -x scripts/proteomicslfq.sh \
  -g \
  -p ~/src/kohlbacher/openms/build/bin/ProteomicsLFQ \
  -r "$results_directory/PXD037527/proteomicslfq" \
  -f "$data_directory/PXD037527/uniprot_reviewed_combined_UP5640_UP2311_2022_08_11.fasta" \
  "${mzml_files[@]}"
