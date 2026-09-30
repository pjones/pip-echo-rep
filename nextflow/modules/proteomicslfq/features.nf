process features {
    input:
    tuple val(id), path(idxml_file)

    output:
    tuple val(id), path("${idxml_file.simpleName}_perc.idXML")

    script:
    """
    PSMFeatureExtractor \
      -in "${idxml_file}" \
      -out "${idxml_file.simpleName}_feat.idXML" \
      -threads ${task.cpus} \
      -debug 0

    PercolatorAdapter \
      -in "${idxml_file.simpleName}_feat.idXML" \
      -out "${idxml_file.simpleName}_perc.idXML" \
      -threads ${task.cpus} \
      -subset_max_train 300000 \
      -decoy_pattern "DECOY_" \
      -post_processing_tdc \
      -score_type pep \
      -debug 0
    """
}
