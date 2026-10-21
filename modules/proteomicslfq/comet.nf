nextflow.enable.types = true

process comet {
    input:
    mzml_file: Path
    fasta_file: Path

    output:
    file("${mzml_file.simpleName}.idXML")

    script:
    def base = mzml_file.simpleName
    def comet_file = "${base}_comet.idXML"
    def fdr_file = "${base}_fdr.idXML"
    def id_filter_file = "${base}.idXML"

    """
    CometAdapter \
      -in "${mzml_file}" \
      -out "${comet_file}" \
      -threads ${task.cpus} \
      -database "${fasta_file}" \
      -instrument low_res \
      -missed_cleavages 2 \
      -min_peptide_length 6 \
      -max_peptide_length 40 \
      -num_hits 1 \
      -num_enzyme_termini fully \
      -enzyme "Trypsin/P" \
      -isotope_error 0/1 \
      -precursor_charge 2:4 \
      -fixed_modifications 'Carbamidomethyl (C)' \
      -variable_modifications 'Acetyl (Protein N-term)' 'Oxidation (M)' \
      -max_variable_mods_in_peptide 3 \
      -precursor_mass_tolerance 20.0 \
      -precursor_error_units ppm \
      -fragment_mass_tolerance 0.5 \
      -fragment_bin_offset 0.4 \
      -minimum_peaks 10 \
      -PeptideIndexing:IL_equivalent \
      -PeptideIndexing:unmatched_action warn \
      -debug 0 \
      -force

    PSMFeatureExtractor \
      -in "${comet_file}" \
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

    FalseDiscoveryRate \
      -in "${base}_percolator.idXML" \
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
