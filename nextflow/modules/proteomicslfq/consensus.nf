process consensus {
    input:
    val files

    output:
    path "${files[0].simpleName}.idXML", emit: idxml

    script:
    def consensus_file = "${files[0].simpleName}_consensus.idXML"
    def fdr_file = "${files[0].simpleName}_fdr.idXML"
    def id_filter_file = "${files[0].simpleName}.idXML"

    """
    ConsensusID \
      -in ${files.join(' ')} \
      -out "${consensus_file}" \
      -per_spectrum \
      -threads ${task.cpus} \
      -algorithm best \
      -filter:min_support 0 \
      -filter:considered_hits 0 \
      -debug 0

    FalseDiscoveryRate \
      -in "${consensus_file}" \
      -out "${fdr_file}" \
      -threads ${task.cpus} \
      -algorithm:add_decoy_peptides \
      -algorithm:add_decoy_proteins \
      -PSM true -protein false

    IDFilter \
      -in "${fdr_file}" \
      -out "${id_filter_file}" \
      -threads ${task.cpus} \
      -score:psm "0.10"
    """
}
