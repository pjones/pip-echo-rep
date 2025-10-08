#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
option_results_dir=
option_fasta_file=

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options] mzML [mzML...]

  -f FILE FASTA file
  -h      This message
  -r DIR  Where to store results

EOF
}

################################################################################
function log() {
  echo "==> " "$@"
}

################################################################################
function bootstrap_metamorpheus() {
  local settings_dir=$1
  local templates_dir=$2
  local calibration_task_file=$3
  local search_task_file=$4

  if [ ! -d "$settings_dir" ]; then
    log "generating MetaMorpheus settings directory"
    MetaMorpheus --mmsettings "$settings_dir" >/dev/null || :

    # Delete all the bullshit:
    find "$settings_dir" \
      -regextype posix-extended \
      -type f \
      \! -regex '.*\.(txt|tsv|xml|mzMl|dat|fasta|mzML|gdb|pdb)(.gz)?$' \
      -delete

    # Don't prompt:
    cat <<EOF >"$settings_dir/settings.toml"
WriteExcelCompatibleTSVs = false
UserHasAgreedToThermoRawFileReaderLicence = true
EOF
  fi

  if [ ! -d "$templates_dir" ]; then
    log "generating default configuration templates for MetaMorpheus"
    mkdir --parents "$templates_dir"

    pushd "$templates_dir" || exit 1
    MetaMorpheus --mmsettings "$settings_dir" -g -o .
    popd || exit 1
  fi

  if [ -n "${calibration_task_file:-}" ] && [ ! -e "$calibration_task_file" ]; then
    cp "$templates_dir/CalibrationTask.toml" "$calibration_task_file"
  fi

  if [ -n "${search_task_file:-}" ] && [ ! -e "$search_task_file" ]; then
    # Modify some of the default settings to match the paper:
    #
    # NOTE: These are the settings that are /different/ than the default
    # settings.
    sed -E \
      -e '/^DoLabelFreeQuantification/s/=.*$/= false/' \
      -e '/^DoLocalizationAnalysis/s/=.*$/= false/' \
      -e '/^WriteIndividualFiles/s/=.*$/= false/' \
      -e '/^QValueThreshold /s/=.*$/= 1.0/' \
      -e '/^PepQValueThreshold /s/=.*$/= 0.01/' \
      -e '/^MinAllowedInternalFragmentLength /s/=.*/= 7/' \
      "$templates_dir/SearchTask.toml" >"$search_task_file"
  fi
}

################################################################################
# Run MetaMorpheus to create PSM ID files.
#
# Arguments:
#
# $1: FASTA file
# $@: Spectra files
function metamorpheus() {
  local fasta_file=$1
  shift

  local spectra_files=()
  local search_task_file=Task1-SearchTask.toml
  local calibration_task_file=Task2-CalibrationTask.toml

  local settings_dir
  settings_dir=$(realpath mmsettings)

  local templates_dir
  templates_dir=$(realpath mmtemplates)

  # Generate settings and configuration:
  bootstrap_metamorpheus \
    "$settings_dir" \
    "$templates_dir" \
    "$calibration_task_file" \
    "$search_task_file"

  # Reuse the search task:
  local task_3_search="${search_task_file//Task1/Task3}"
  cp "$search_task_file" "$task_3_search"

  # Do the identification:
  if [ ! -e Task1SearchTask/results.txt ]; then
    log "running MetaMorpheus"
    time MetaMorpheus \
      --mmsettings "$settings_dir" \
      -t "$search_task_file" "$calibration_task_file" "$task_3_search" \
      -d "$fasta_file" \
      -s "${spectra_files[@]}" \
      -o . |
      tee "metamorpheus.log"
  fi
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
  while getopts "f:hr:" o; do
    case "${o}" in
    f)
      option_fasta_file=$OPTARG
      ;;

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

  metamorpheus "$option_fasta_file" "$@"
  flashlfq Task3SearchTask/AllPSMs.psmtsv Task2CalibrationTask
}

################################################################################
main "$@"
