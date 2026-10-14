process msgfplus {
    input:
    path mzml_file
    path fasta_file

    output:
    tuple val(mzml_file.simpleName), path("${mzml_file.simpleName}_msgfplus.idXML")

    script:
    """
    msgfplus=\$(which msgfplus)
    jar=\$(realpath "\$(dirname "\$msgfplus")/../share/msgfplus/MSGFPlus.jar")

    MSGFPlusAdapter \
      -protocol automatic \
      -in "${mzml_file}" \
      -out "${mzml_file.simpleName}_msgfplus.idXML" \
      -executable "\$jar" \
      -threads ${task.cpus} \
      -java_memory 6144 \
      -database "${fasta_file}" \
      -instrument low_res \
      -matches_per_spec 1 \
      -min_precursor_charge 2 \
      -max_precursor_charge 4 \
      -min_peptide_length 6 \
      -max_peptide_length 40 \
      -max_missed_cleavages 2 \
      -isotope_error_range 0,1 \
      -enzyme "Trypsin/P" \
      -tryptic fully \
      -precursor_mass_tolerance 20.0 \
      -precursor_error_units ppm \
      -fixed_modifications 'Carbamidomethyl (C)' \
      -variable_modifications 'Acetyl (Protein N-term)' 'Oxidation (M)' \
      -max_mods 3 \
      -PeptideIndexing:IL_equivalent \
      -PeptideIndexing:unmatched_action warn \
      -debug 0
    """
}
