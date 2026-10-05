process flashlfq {
    input:
    tuple path(psmtsv), path(peptides_tsv), path(calibration_task_dir)

    output:
    tuple path("QuantifiedPeaks.tsv"),
          path("QuantifiedPeptides.tsv"),
          path("QuantifiedProteins.tsv"),
          path("flashlfq.log")

    script:
    """
    flashlfq.sh -d ${calibration_task_dir} -p "${psmtsv}"
    """
}
