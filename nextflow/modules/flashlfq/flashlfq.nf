nextflow.enable.types = true

process flashlfq {
    input:
    record(
      t3_psms: Path,
      t2_cdir: Path,
    )

    output:
    record(
      peaks: file("QuantifiedPeaks.tsv"),
      peptides: file("QuantifiedPeptides.tsv"),
      proteins: file("QuantifiedProteins.tsv")
    )

    script:
    """
    flashlfq.sh -d ${t2_cdir} -p "${t3_psms}"
    """
}
