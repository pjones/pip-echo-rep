nextflow.enable.types = true

process comet {
    input:
    mzml_file: Path
    fasta_file: Path

    output:
    record(
        base: mzml_file.simpleName,
        idxml: file("${mzml_file.simpleName}_comet.idXML")
    )

    script:
    """
    CometAdapter \
      -in "${mzml_file}" \
      -out "${mzml_file.simpleName}_comet.idXML" \
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
    """
}
