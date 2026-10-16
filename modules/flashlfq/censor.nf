nextflow.enable.types = true

process censor {
    input:
    mzml_files: List<Path>
    record(
        t1_psms: Path,
        t1_peps: Path,
    )

    output:
    record(
        mzml: files("*-censored.mzML"),
        psms: file("censored.psmtsv")
    )

    script:
    """
    censor.py --tool flashlfq ${t1_peps} ${t1_psms} ${mzml_files.join(' ')}
    """
}
