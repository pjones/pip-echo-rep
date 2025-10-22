#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
option_results_dir=
option_fasta_file=
option_generate_decoys=1
option_decoy_prefix="DECOY_"
option_design_file=
option_proteomicslfq=ProteomicsLFQ
option_use_pip_echo=0
option_use_gdb=0
option_gdb_flags=()

################################################################################
declare -a id_xml_files

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options] mzML [mzML...]

  -d FILE Design file
  -D      Turn off decoy generation
  -e      Use PIP-ECHO
  -f FILE FASTA file
  -g      Enable debugging via GDB
  -G OPT  Pass OPT to GDB
  -h      This message
  -p BIN  The ProteomicsLFQ binary to use
  -r DIR  Where to store results

EOF
}

################################################################################
# Given the name of an ID file, generate the Percolator output name
# for that file.
function make_perc_file() {
  local id_file=$1
  echo "$(dirname "$id_file")/$(basename "$id_file" .idXML)_perc.idXML"
}

################################################################################
function id_with_msgf() {
  # This is a global variable:
  id_with_msgf_res=()

  # OpenMS runs the JAR directly :(
  local jar
  local msgfplus
  msgfplus=$(which msgfplus)
  jar=$(realpath "$(dirname "$msgfplus")/../share/msgfplus/MSGFPlus.jar")

  local index_dir="msgf_index"
  mkdir -p "$index_dir"

  local msgf_index
  msgf_index="$index_dir/$(basename "$option_fasta_file" .fasta).cseq"

  local fasta_file=$option_fasta_file
  local buildsa_flags=()

  if [ "$option_generate_decoys" -eq 1 ]; then
    buildsa_flags+=("-decoy" "$option_decoy_prefix")
    buildsa_flags+=("-tda" "1")
    msgf_index="$index_dir/$(basename "$option_fasta_file" .fasta).revCat.cseq"
    option_fasta_file=$(realpath "$index_dir/$(basename "$option_fasta_file" .fasta).revCat.fasta")
  else
    buildsa_flags+=("-tda" "0")
  fi

  if [ ! -e "$msgf_index" ]; then
    # https://msgfplus.github.io/msgfplus/BuildSA.html
    java \
      -Xms6144m \
      -Xmx6144m \
      -cp "$jar" \
      edu.ucsd.msjava.msdbsearch.BuildSA \
      -d "$fasta_file" \
      -o "$index_dir" \
      -debug 0 \
      "${buildsa_flags[@]}"
  fi

  local msgf_id_dir="msgf_ids"
  mkdir -p "$msgf_id_dir"

  for mzml_file in "$@"; do
    local id_file
    id_file="$msgf_id_dir/$(basename "$mzml_file" .mzML).idXML"
    id_with_msgf_res+=("$id_file")

    if [ ! -e "$id_file" ]; then
      # FIXME: Need to pass correct PSM and enzyme information.
      MSGFPlusAdapter \
        -protocol automatic \
        -in "$mzml_file" \
        -out "$id_file" \
        -executable "$jar" \
        -threads 2 \
        -java_memory 6144 \
        -database "$option_fasta_file" \
        -instrument low_res \
        -matches_per_spec 1 \
        -min_precursor_charge 2 \
        -max_precursor_charge 4 \
        -min_peptide_length 6 \
        -max_peptide_length 40 \
        -max_missed_cleavages 2 \
        -isotope_error_range 0,1 \
        -enzyme "Trypsin/P" \
        -tryptic fully \
        -precursor_mass_tolerance 20.0 \
        -precursor_error_units ppm \
        -fixed_modifications 'Carbamidomethyl (C)' \
        -variable_modifications 'Acetyl (Protein N-term)' 'Oxidation (M)' \
        -max_mods 3 \
        -PeptideIndexing:IL_equivalent \
        -PeptideIndexing:unmatched_action warn \
        -debug 0
    fi
  done
}

################################################################################
function id_with_comet() {
  # Global variable:
  id_with_comet_res=()

  local comet_id_dir="comet_ids"
  mkdir -p "$comet_id_dir"

  for mzml_file in "$@"; do
    local id_file
    id_file="$comet_id_dir/$(basename "$mzml_file" .mzML).idXML"
    id_with_comet_res+=("$id_file")

    if [ ! -e "$id_file" ]; then
      CometAdapter \
        -in "$mzml_file" \
        -out "$id_file" \
        -threads 2 \
        -database "$option_fasta_file" \
        -instrument low_res \
        -missed_cleavages 2 \
        -min_peptide_length 6 \
        -max_peptide_length 40 \
        -num_hits 1 \
        -num_enzyme_termini fully \
        -enzyme "Trypsin/P" \
        -isotope_error 0/1 \
        -precursor_charge 2:4 \
        -fixed_modifications 'Carbamidomethyl (C)' \
        -variable_modifications 'Acetyl (Protein N-term)' 'Oxidation (M)' \
        -max_variable_mods_in_peptide 3 \
        -precursor_mass_tolerance 20.0 \
        -precursor_error_units ppm \
        -fragment_mass_tolerance 0.5 \
        -fragment_bin_offset 0.4 \
        -minimum_peaks 10 \
        -PeptideIndexing:IL_equivalent \
        -PeptideIndexing:unmatched_action warn \
        -debug 0 \
        -force
    fi
  done
}

################################################################################
function postprocess_id_files() {
  for id_file in "$@"; do
    local feature_file
    feature_file="$(dirname "$id_file")/$(basename "$id_file" .idXML)_feat.idXML"

    if [ ! -e "$feature_file" ]; then
      PSMFeatureExtractor \
        -in "$id_file" \
        -out "$feature_file" \
        -threads 2 \
        -debug 0
    fi

    local perc_file
    perc_file="$(make_perc_file "$id_file")"

    if [ ! -e "$perc_file" ]; then
      PercolatorAdapter \
        -in "$feature_file" \
        -out "$perc_file" \
        -threads 2 \
        -subset_max_train 300000 \
        -decoy_pattern "$option_decoy_prefix" \
        -post_processing_tdc \
        -score_type pep \
        -debug 0
    fi
  done
}

################################################################################
function join_id_files() {
  # Global variable:
  join_id_files_res=()

  local msgf_perc_file
  local comet_perc_file
  local consensus_file

  local consensus_dir="consensus_id"
  mkdir -p "$consensus_dir"

  for i in "${!id_with_msgf_res[@]}"; do
    consensus_file=$consensus_dir/$(
      basename "${id_with_msgf_res[i]}"
    )

    join_id_files_res+=("$consensus_file")
    msgf_perc_file=$(make_perc_file "${id_with_msgf_res[i]}")
    comet_perc_file=$(make_perc_file "${id_with_comet_res[i]}")

    if [ ! -e "$consensus_file" ]; then
      ConsensusID \
        -in "$msgf_perc_file" "$comet_perc_file" \
        -out "$consensus_file" \
        -per_spectrum \
        -threads 1 \
        -algorithm best \
        -filter:min_support 0 \
        -filter:considered_hits 0 \
        -debug 0
    fi
  done
}

################################################################################
function filter_consensus_id_files() {
  # Global variable:
  filter_consensus_id_files_res=()

  local fdr_file
  local id_filter_file

  for consensus_file in "$@"; do
    fdr_file=$(dirname "$consensus_file")/$(basename "$consensus_file" .idXML)_fdr.idXML
    id_filter_file=$(dirname "$consensus_file")/$(basename "$consensus_file" .idXML)_filter.idXML
    filter_consensus_id_files_res+=("$(realpath "$id_filter_file")")

    if [ ! -e "$fdr_file" ]; then
      FalseDiscoveryRate \
        -in "$consensus_file" \
        -out "$fdr_file" \
        -threads 2 \
        -algorithm:add_decoy_peptides \
        -algorithm:add_decoy_proteins \
        -PSM true -protein false
    fi

    if [ ! -e "$id_filter_file" ]; then
      IDFilter \
        -in "$fdr_file" \
        -out "$id_filter_file" \
        -threads 2 \
        -score:psm "0.10"
    fi
  done
}

################################################################################
function identify_peptides() {
  id_with_msgf "$@"
  postprocess_id_files "${id_with_msgf_res[@]}"

  id_with_comet "$@"
  postprocess_id_files "${id_with_comet_res[@]}"

  join_id_files
  filter_consensus_id_files "${join_id_files_res[@]}"
  id_xml_files=("${filter_consensus_id_files_res[@]}")

  # Clean up some memory
  filter_consensus_id_files_res=()
  join_id_files_res=()
  id_with_msgf_res=()
  id_with_comet_res=()
}

################################################################################
function quantify_peptides() {
  local mz_tab_file="proteomicslfq.mzTab"
  local proteomicslfq_flags=()
  local command=$option_proteomicslfq
  local options=()

  if [ -n "$option_design_file" ]; then
    proteomicslfq_flags+=("-design" "$option_design_file")
  fi

  if [ "$option_use_pip_echo" -eq 1 ]; then
    proteomicslfq_flags+=("-pip_echo" "true")
  fi

  if [ "$option_use_gdb" -eq 1 ]; then
    command=gdb
    options=("${option_gdb_flags[@]}")
    options+=("--args" "$option_proteomicslfq")
  fi

  if [ "$option_use_gdb" -eq 1 ] || [ ! -e "$mz_tab_file" ]; then
    "$command" "${options[@]}" \
      -threads 2 \
      -in "$@" \
      -ids "${id_xml_files[@]}" \
      -fasta "$option_fasta_file" \
      -protein_inference aggregation \
      -quantification_method feature_intensity \
      -targeted_only false \
      -feature_with_id_min_score 0.10 \
      -mass_recalibration false \
      -Seeding:intThreshold 1000 \
      -protein_quantification unique_peptides \
      -alignment_order star \
      -PeptideQuantification:quantify_decoys \
      -psmFDR 0.01 \
      -proteinFDR 0.01 \
      -picked_proteinFDR true \
      -out "$mz_tab_file" \
      -out_cxml proteomicslfq.consensusXML \
      -out_msstats msstats.csv \
      -out_triqler triqler.tsv \
      -debug 0 \
      "${proteomicslfq_flags[@]}"
  fi
}

################################################################################
function main() {
  while getopts "d:Def:gG:hp:r:" o; do
    case "${o}" in
    d)
      option_design_file=$OPTARG

      if [ ! -e "$option_design_file" ]; then
        echo >&2 "ERROR: design file doesn't exist: $option_design_file"
        exit 1
      fi
      ;;

    D)
      option_generate_decoys=0
      ;;

    e)
      option_use_pip_echo=1
      ;;

    f)
      option_fasta_file=$OPTARG
      ;;

    g)
      option_use_gdb=1
      ;;

    G)
      option_gdb_flags+=("$OPTARG")
      ;;

    h)
      usage
      exit
      ;;

    p)
      option_proteomicslfq=$OPTARG
      ;;

    r)
      option_results_dir=$OPTARG
      ;;

    *)
      exit 1
      ;;
    esac
  done

  shift $((OPTIND - 1))

  if [ -z "${option_fasta_file:-}" ] ||
    [ ! -e "$option_fasta_file" ]; then
    echo >&2 "ERROR: FASTA file not given or does not exist: ${option_fasta_file:-}"
    exit 1
  fi

  if [ -z "${option_results_dir:-}" ]; then
    echo >&2 "ERROR: no results directory option given"
    exit 1
  fi

  if [ $# -lt 1 ]; then
    echo >&2 "ERROR: missing spectrum files"
    exit 1
  fi

  mkdir --parents "$option_results_dir"
  pushd "$option_results_dir" || exit 1

  identify_peptides "$@"
  quantify_peptides "$@"
}

################################################################################
main "$@"
