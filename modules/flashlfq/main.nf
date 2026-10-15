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
include { analysis } from "../analysis.nf"

workflow flashlfq {
    take:
    mbr_fdr: Float
    mzml_files: List<Path>
    fasta_files: List<Path>

    main:
    id_before = mm_before(mzml_files, fasta_files)
    quant_before = lfq_before(mbr_fdr, id_before)

    censored = censor(mzml_files, id_before)
    censored_mzml = censored.map { c -> c.mzml.toSorted() }

    id_after = mm_after(censored_mzml, fasta_files)
    quant_after = lfq_after(mbr_fdr, id_after)

    analyzed = analysis(
        "flashlfq",
        params.data_dir,
        censored.map { r -> r.psms },
        quant_before.map { r -> r.peaks.parent },
        quant_after.map { r -> r.peaks.parent },
        mbr_fdr
    )

    emit:
    censored: Value<Record> = censored
    id_before: Value<Record> = id_before
    quant_before: Value<Record> = quant_before
    id_after: Value<Record> = id_after
    quant_after: Value<Record> = quant_after
    analyzed: Channel<Path> = analyzed
}
