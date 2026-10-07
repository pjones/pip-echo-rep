nextflow.enable.types = true

process censor {
    input:
    mzml_files: List<Path>
    record(
        t1_psms: Path,
        t1_peps: Path,
    )

    output:
    files("*-censored.mzML")

    script:
    """
    censor.py --peptides ${t1_peps} --psms ${t1_psms} --mzml ${mzml_files.join(' ')}
    """
}
