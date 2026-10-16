nextflow.enable.types = true

process censor {
    input:
    files: List<Path>

    output:
    record(
        mzml: files("*-censored.mzML"),
        censored: file("censored.tsv"),
    )

    script:
    """
    censor.py --tool proteomicslfq .
    """
}
