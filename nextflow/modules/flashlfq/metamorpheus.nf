process metamorpheus {
    input:
    path mzml_files
    path fast_files

    output:
    path "Task3SearchTask/AllPSMs.psmtsv", emit: "psmtsv"
    path "Task2CalibrationTask", emit: "calibration_task_dir"

    script:
    """
    metamorpheus.sh -f ${fast_files.join(' -f ')} ${mzml_files.join(' ')}
    """
}
