process toMzML {
    input:
    path raw_file

    output:
    path "${raw_file.simpleName}.mzML"

    script:
    """
    ThermoRawFileParser \
      --input="${raw_file}" \
      --output="${raw_file.simpleName}.mzML" \
      --format=2 \
      --msLevel=1-
    """
}
