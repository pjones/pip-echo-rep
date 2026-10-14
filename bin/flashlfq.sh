#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
option_ppm_tolerance="10.0"
option_isotope_ppm_tolerance="5.0"
option_num_isotopes_required="2"
option_mbr_rt_window="1.5" # Minutes
option_fdr="0.01"

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

  -d DIR  Directory containing Metamorpheus output mzML files
  -h      This message
  -p FILE The PSMs TSV file

  -F FLOAT MBR FDR
  -I FLOAT Isotope PPM tolerance
  -N INT   Number of isotopes required
  -P FLOAT PPM Tolerance
  -W FLOAT MBR RT Window

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
      --rea \
      --ath \
      --out . \
      --idt "$fixed_psms_file" \
      --rep "$mzml_dir" \
      --ppm "${option_ppm_tolerance}" \
      --iso "${option_isotope_ppm_tolerance}" \
      --nis "${option_num_isotopes_required}" \
      --mbr \
      --mrt "${option_mbr_rt_window}" \
      --pipfdr "${option_fdr}" |
      tee "flashlfq.log"
  fi
}

################################################################################
function main() {
  local psms_file
  local mzml_dir

  while getopts "d:hp:F:P:I:N:W:" o; do
    case "${o}" in
    d)
      mzml_dir=$OPTARG
      ;;

    h)
      usage
      exit
      ;;

    F)
      option_fdr=$OPTARG
      ;;

    I)
      option_isotope_ppm_tolerance=$OPTARG
      ;;

    N)
      option_num_isotopes_required=$OPTARG
      ;;

    p)
      psms_file=$OPTARG
      ;;

    P)
      option_ppm_tolerance=$OPTARG
      ;;

    W)
      option_mbr_rt_window=$OPTARG
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
