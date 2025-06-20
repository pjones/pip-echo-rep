#!/usr/bin/env bash
data_directory=data

function download_pride_ftp_url() {
  local url=$1
  local id
  local cache_file

  id=$(basename "$url")
  cache_file="$data_directory/$id/.ls_cache"

  # Ensure there's a slash on the end of the URL:
  if [ "${url: -1:1}" != "/" ]; then
    url="$url/"
  fi

  mkdir --parents "$data_directory/$id"

  # Download and cache the directory listing:
  if [ ! -e "$cache_file" ]; then
    curl --list-only "$url" >"$cache_file"
  fi

  # Make sure each of the files has been downloaded:
  while read -r file; do
    path="$data_directory/$id/$file"

    if [ ! -e "$path" ]; then
      curl --output "$path" "${url}${file}"
    fi
  done <"$cache_file"

  # Report the size of the directory
  du --human-readable --summarize "$data_directory/$id"
}

main() {
  if [ $# -gt 0 ]; then
    data_directory=$1
  fi

  download_pride_ftp_url "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD014415/"
  download_pride_ftp_url "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/05/PXD003881/"
  download_pride_ftp_url "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/07/PXD037527/"
}

main "$@"
