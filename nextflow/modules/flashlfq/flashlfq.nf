process flashlfq {
    input:
    path psmtsv
    path calibration_task_dir

    output:
    path "flashlfq_AllPSMs.psmtsv", emit: "psmtsv"
    path "flashlfq.log", emit: "log"

    script:
    """
    flashlfq.sh -d ${calibration_task_dir} -p "${psmtsv}"
    """
}
