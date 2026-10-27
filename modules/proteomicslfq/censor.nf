nextflow.enable.types = true

process censor {
    input:
    files: List<Path>

    output:
    record(
        mzml: files("*-censored.mzML"),
        censored: file("censored.tsv"),
        log: file("censor.log"),
    )

    script:
    """
    censor.py \
        --log censor.log \
        --tool proteomicslfq \
        .
    """
}
