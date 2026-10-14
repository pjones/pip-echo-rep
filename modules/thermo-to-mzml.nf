nextflow.enable.types = true

process toMzML {
    input:
    raw_file: Path

    output:
    file("${raw_file.simpleName}.mzML")

    script:
    """
    ThermoRawFileParser \
      --input="${raw_file}" \
      --output="${raw_file.simpleName}.mzML" \
      --format=2 \
      --msLevel=1-
    """
}
