process proteomicslfq {
    input:
    val mzml_files
    val idxml_files
    path fasta_file

    output:
    tuple path("proteomicslfq.mzTab"),
          path("proteomicslfq.consensusXML"),
          path("msstats.csv")

    script:
    """
    ProteomicsLFQ \
      -in ${mzml_files.join(' ')} \
      -ids ${idxml_files.join(' ')} \
      -fasta "${fasta_file}" \
      -pip_echo true \
      -threads ${task.cpus} \
      -protein_inference aggregation \
      -quantification_method feature_intensity \
      -targeted_only false \
      -feature_with_id_min_score 0.10 \
      -mass_recalibration false \
      -Seeding:intThreshold 1000 \
      -protein_quantification unique_peptides \
      -alignment_order star \
      -PeptideQuantification:quantify_decoys \
      -psmFDR 0.01 \
      -proteinFDR 0.01 \
      -picked_proteinFDR true \
      -out "proteomicslfq.mzTab" \
      -out_cxml proteomicslfq.consensusXML \
      -out_msstats msstats.csv \
      -debug 0
    """
}
