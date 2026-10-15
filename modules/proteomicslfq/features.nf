nextflow.enable.types = true

process features {
    input:
    record(
        base: String,
        idxml: Path,
    )

    output:
    record(
        base: base,
        feature: file("${base}_feature.idXML"),
        percolator: file("${base}_percolator.idXML"),
    )

    script:
    """
    PSMFeatureExtractor \
      -in "${idxml}" \
      -out "${base}_feature.idXML" \
      -threads ${task.cpus} \
      -debug 0

    PercolatorAdapter \
      -in "${base}_feature.idXML" \
      -out "${base}_percolator.idXML" \
      -threads ${task.cpus} \
      -subset_max_train 300000 \
      -decoy_pattern "DECOY_" \
      -post_processing_tdc \
      -score_type pep \
      -debug 0
    """
}
