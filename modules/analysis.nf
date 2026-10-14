nextflow.enable.types = true

process analysis {
    input:
    tool_name: String
    data_dir: Path
    censored_manifest: Path
    orig_path: Path
    censored_path: Path
    mbr_fdr: Float

    output:
    record(
        json: file("analysis-${mbr_fdr}.json"),
        log: file("analysis-${mbr_fdr}.log"),
    )

    script:
    """
    analysis.py \
        --output "analysis-${mbr_fdr}.json" \
        --log "analysis-${mbr_fdr}.log" \
        --tool "${tool_name}" \
        --id "${params.dataset}" \
        --qval "${mbr_fdr}" \
        --human-file-pattern "${params.human_file_pattern}" \
        --foreign-species "${params.foreign_species}" \
        --rt-delta "${params.rt_delta}" \
        --data "${data_dir}" \
        --censored-manifest "${censored_manifest}" \
        "${orig_path}" "${censored_path}"
    """
}
