nextflow.enable.types = true

include {
    metamorpheus as mm_before ;
    metamorpheus as mm_after
} from "./metamorpheus.nf"

include {
    flashlfq as lfq_before ;
    flashlfq as lfq_after
} from "./flashlfq.nf"

include { censor } from "./censor.nf"

workflow flashlfq {
    take:
    mzml_files: List<Path>
    fasta_files: List<Path>

    main:
    id_before = mm_before(mzml_files, fasta_files)
    quant_before = lfq_before(id_before)

    censored_mzml = censor(mzml_files, id_before).map { c -> c.toSorted() }

    id_after = mm_after(censored_mzml, fasta_files)
    quant_after = lfq_after(id_after)

    emit:
    censored_mzml: Value<List<Path>> = censored_mzml
    quant_before: Value<Record> = quant_before
    quant_after: Value<Record> = quant_after
}
