include { flashlfq } from "./flashlfq/main.nf"
include { proteomicslfq } from "./proteomicslfq/main.nf"

workflow analysis {
    take:
    mzml_files
    fasta_files

    main:
    flashlfq(mzml_files, fasta_files)
    proteomicslfq(mzml_files, fasta_files)

    emit:
    flashlfq_id = flashlfq.out.id
    flashlfq_quant = flashlfq.out.quant

    proteomicslfq = proteomicslfq.out
}
