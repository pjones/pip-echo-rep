nextflow.enable.types = true

process analysis {
    input:
    tool_name: String
    data_dir: Path
    censored_manifest: Path
    record(
        orig_path: Path,
        censored_path: Path,
        mbr_fdr: Float,
    )

    output:
    record(
        mbr_fdr: mbr_fdr,
        json: file("analysis.json"),
        log: file("analysis.log"),
        csv: files("*.csv"),
    )

    script:
    """
    analysis.py \
        --output "analysis.json" \
        --log "analysis.log" \
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
