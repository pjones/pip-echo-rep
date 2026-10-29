nextflow.enable.types = true

process proteomicslfq {
    input:
    mbr_fdr: Float
    fasta_file: Path
    mzml_files: List<Path>
    idxml_files: List<Path>

    output:
    record(
        mbr_fdr: mbr_fdr,
        mztab: file("proteomicslfq.mzTab"),
        consensus: file("proteomicslfq.consensusXML"),
        msstats: file("msstats.csv"),
        trafoXML: files("*.trafoXML"),
        log: file("proteomicslfq.log"),
    )

    script:
    def pip_echo_enabled = mbr_fdr < 1.0 ? "true" : "false"

    """
    ProteomicsLFQ \
      -in ${mzml_files.join(' ')} \
      -ids ${idxml_files.join(' ')} \
      -fasta "${fasta_file}" \
      -pip_echo ${pip_echo_enabled} \
      -PipEcho:fdr "${mbr_fdr}" \
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
      -out proteomicslfq.mzTab \
      -out_cxml proteomicslfq.consensusXML \
      -out_msstats msstats.csv \
      -debug 1000 |
      tee "proteomicslfq.log"
    """
}
