#!/usr/bin/env bash

################################################################################
data_directory=$(realpath "data")
results_directory=$(realpath "results")/$(basename "$0" .sh)
option_force=0

################################################################################
usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

  -h      This message
  -f      Force run

EOF
}

################################################################################
function main() {
  # Option arguments are in $OPTARG
  while getopts "fh" o; do
    case "${o}" in
    f)
      option_force=1
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

  if [ -e "$results_directory" ] && [ "$option_force" -eq 0 ]; then
    echo >&2 "ERROR: results directory already exists, use -f"
    return 1
  fi

  mkdir --parents "$results_directory"
  pushd "$results_directory" || exit 1

  # We need the ID files from the following ZIP files:
  local metamorpheus=SingleCell_SearchResults-MetaMorpheus

  if [ ! -d "$metamorpheus" ]; then
    unzip "$data_directory/PXD057758/SingleCell_SearchResults-MetaMorpheus.zip"
  fi

  if [ ! -e LicenseAgreements.toml ]; then
    echo "HasAcceptedThermoLicence = true" >LicenseAgreements.toml
  fi

  # FIXME: what command line argument sets this?
  # DonorQValueThreshold = 0.002

  set -x

  time FlashLFQ \
    --rea --ath \
    --out "$(realpath .)" \
    --idt ./SingleCell_SearchResults-MetaMorpheus/Task1-SearchTask-Censored/AllPSMs.psmtsv \
    --rep "$data_directory/PXD037527" \
    --ppm 10.0 --iso 5.0 --nis 2 \
    --mbr --mrt 1.5 --pipfdr 0.01 |
    tee "$(date +%s).log"
}

################################################################################
main "$@"
