#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

  -d DIR  Directory containing Metamorpheus output mzML files
  -h      This message
  -p FILE The PSMs TSV file

EOF
}

################################################################################
function log() {
  echo "==> " "$@"
}

################################################################################
function flashlfq() {
  local psms_file=$1
  local mzml_dir=$2

  # FlashLFQ and MetaMorpheus don't agree on the column names anymore :(
  local fixed_psms_file
  fixed_psms_file=flashlfq_$(basename "$psms_file")

  if [ ! -e "$fixed_psms_file" ]; then
    sed \
      -e '1s/\tSpectrum Match Count\t/\tPSM Count (unambiguous, <0.01 q-value)\t/' \
      -e '1s/\tMonoisotopic Mass\t/\tPeptide Monoisotopic Mass\t/' \
      -e '1s/\tAccession\t/\tProtein Accession\t/' \
      -e '1s/\tName\t/\tProtein Name\t/' \
      -e '1s/\tDescription\t/\tPeptide Description\t/' \
      -e '1s/\tStart and End Residues In Full Sequence\t/\tStart and End Residues In Protein\t/' \
      -e '1s/\tPrevious Residue\t/\tPrevious Amino Acid\t/' \
      -e '1s/\tNext Residue\t/\tNext Amino Acid\t/' \
      "$psms_file" >"$fixed_psms_file"
  fi

  if [ ! -e QuantifiedPeptides.tsv ]; then
    log "running FlashLFQ"
    time FlashLFQ \
      --rea --ath \
      --out . --idt "$fixed_psms_file" --rep "$mzml_dir" \
      --ppm 10.0 --iso 5.0 --nis 2 \
      --mbr --mrt 1.5 --pipfdr 0.01 |
      tee "flashlfq.log"
  fi
}

################################################################################
function main() {
  local psms_file
  local mzml_dir

  while getopts "d:hp:" o; do
    case "${o}" in
    d)
      mzml_dir=$OPTARG
      ;;

    h)
      usage
      exit
      ;;

    p)
      psms_file=$OPTARG
      ;;

    *)
      exit 1
      ;;
    esac
  done

  shift $((OPTIND - 1))

  if [ $# -gt 0 ]; then
    echo >&2 "ERROR: unexpected arguments: $*"
    exit 1
  fi

  if [ -z "${psms_file:-}" ] || [ ! -e "$psms_file" ]; then
    echo >&2 "ERROR: missing PSMs file"
    exit 1
  fi

  if [ -z "${mzml_dir:-}" ] || [ ! -e "$mzml_dir" ]; then
    echo >&2 "ERROR: missing Metamorpheus mzML directory"
    exit 1
  fi

  flashlfq "$psms_file" "$mzml_dir"
}

################################################################################
main "$@"
