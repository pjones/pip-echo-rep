#!/usr/bin/env bash

################################################################################
set -eu
set -o pipefail

################################################################################
option_fasta_files=()

################################################################################
function usage() {
  cat <<EOF
Usage: $(basename "$0") [options] mzML [mzML]

  -h       This message
  -f FILE  Add a FASTA file

EOF
}

################################################################################
function log() {
  echo "==> " "$@"
}

################################################################################
function bootstrap() {
  local settings_dir=$1
  local templates_dir=$2
  local calibration_task_file=$3
  local search_task_file=$4

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

  log "generating default configuration templates for MetaMorpheus"
  mkdir --parents "$templates_dir"

  pushd "$templates_dir" || exit 1
  MetaMorpheus --mmsettings "$settings_dir" -g -o .
  popd || exit 1

  cp "$templates_dir/CalibrationTask.toml" "$calibration_task_file"

  # Modify some of the default settings to match the paper:
  #
  # NOTE: These are the settings that are /different/ than the default
  # settings.
  tomlq --toml-output '
      .SearchParameters.DoLocalizationAnalysis |= false |
      .SearchParameters.DoLabelFreeQuantification |= false |
      .SearchParameters.WriteIndividualFiles |= false
    ' <"$templates_dir/SearchTask.toml" >"$search_task_file"
}

################################################################################
# Run MetaMorpheus to create PSM ID files.
#
# Arguments:
#
# $@: Spectra files
function metamorpheus() {
  local spectra_files=("$@")
  local search_task_file=Task1-SearchTask.toml
  local calibration_task_file=Task2-CalibrationTask.toml

  local settings_dir
  settings_dir=$(realpath mmsettings)

  local templates_dir
  templates_dir=$(realpath mmtemplates)

  # Generate settings and configuration:
  bootstrap \
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
      -d "${option_fasta_files[@]}" \
      -s "${spectra_files[@]}" \
      -o . |
      tee "metamorpheus.log"
  fi
}

################################################################################
function main() {
  while getopts "f:h" o; do
    case "${o}" in
    f)
      option_fasta_files+=("$OPTARG")
      ;;

    h)
      usage
      exit
      ;;

    *)
      exit 1
      ;;
    esac
  done

  shift $((OPTIND - 1))

  if [ $# -eq 0 ]; then
    echo >&2 "ERROR: missing mzML files"
    exit 1
  fi

  if [ "${#option_fasta_files[@]}" -eq 0 ]; then
    echo >&2 "ERROR: no FASTA files given"
    exit 1
  fi

  metamorpheus "$@"
}

################################################################################
main "$@"
