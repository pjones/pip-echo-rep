nextflow.enable.types = true

process censor {
    input:
    mzml_files: List<Path>
    record(
        t3_psms: Path,
        t3_peps: Path,
    )

    output:
    record(
        mzml: files("*-censored.mzML"),
        psms: file("censored.tsv"),
        log: file("censor.log"),
    )

    script:
    """
    censor.py \
        --log censor.log \
        --tool flashlfq \
        ${t3_psms} ${mzml_files.join(' ')}
    """
}
