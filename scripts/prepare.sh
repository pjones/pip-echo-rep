#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
option_results_dir=
option_mzml_index="mzml_files.lst"

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options] raw [raw...]

  -h      This message
  -r DIR  Where to store results

EOF
}

################################################################################
# Convert spectra files to the appropriate format.
function convert_spectra_files() {
  local spectra_files=()

  for file in "$@"; do
    case "${file##*.}" in
    raw)
      # Convert raw to mzML
      output=$(basename "$file" .raw).mzML

      if [ ! -e "$output" ]; then
        ThermoRawFileParser \
          --input="$file" \
          --output="$output" \
          --format=2 \
          --msLevel=1-
      fi

      spectra_files+=("$output")
      ;;

    *)
      # Just link the file here:
      output=$(basename "$file")
      spectra_files+=("$output")

      ln \
        --symbolic \
        --relative \
        --force \
        --no-dereference \
        --logical \
        "$file" "$output"
      ;;
    esac
  done

  rm --force "$option_results_dir/$option_mzml_index"

  for file in "${spectra_files[@]}"; do
    realpath "$file" >>"$option_results_dir/$option_mzml_index"
  done
}

################################################################################
function main() {
  while getopts "hr:" o; do
    case "${o}" in
    h)
      usage
      exit
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

  if [ $# -lt 1 ]; then
    echo >&2 "ERROR: missing spectrum files"
    exit 1
  fi

  mkdir --parents "$option_results_dir"
  pushd "$option_results_dir" || exit 1

  convert_spectra_files "$@"
}

################################################################################
main "$@"
