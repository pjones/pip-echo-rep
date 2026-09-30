include { metamorpheus } from "./metamorpheus.nf"
include { flashlfq as lfq } from "./flashlfq.nf"

workflow flashlfq {
    take:
    mzml_files
    fasta_files

    main:
    metamorpheus(
        mzml_files,
        fasta_files,
    )

    lfq(
        metamorpheus.out.psmtsv,
        metamorpheus.out.calibration_task_dir,
    )
}
