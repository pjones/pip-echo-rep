nextflow.enable.types = true

process metamorpheus {
    input:
    mzml_files: List<Path>
    fast_files: List<Path>

    output:
    record(
      t1_psms: file("Task1SearchTask/AllPSMs.psmtsv"),
      t1_peps: file("Task1SearchTask/AllPeptides.psmtsv"),
      t2_cdir: file("Task2CalibrationTask"),
      t3_psms: file("Task3SearchTask/AllPSMs.psmtsv"),
      t3_peps: file("Task3SearchTask/AllPeptides.psmtsv"),
    )

    script:
    """
    metamorpheus.sh -f ${fast_files.join(' -f ')} ${mzml_files.join(' ')}
    """
}
