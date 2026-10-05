process metamorpheus {
    input:
    path mzml_files
    path fast_files

    output:
    tuple path("Task3SearchTask/AllPSMs.psmtsv"),
          path("Task3SearchTask/AllPeptides.psmtsv"),
          path("Task2CalibrationTask")

    script:
    """
    metamorpheus.sh -f ${fast_files.join(' -f ')} ${mzml_files.join(' ')}
    """
}
