nextflow.enable.types = true

include { comet } from "./comet.nf"
include { proteomicslfq as lfq } from "./proteomicslfq.nf"

// Emit the files stored in channel as a sorted list.
def as_sorted_files(ch) {
    ch.collect().map { files -> files.toSorted() {f -> f.simpleName} }
}

workflow id_quant {
    take:
    mbr_fdr: Channel<Float>
    fasta: Value<Path>
    mzml_files: Channel<Path>

    main:
    id_files = comet(mzml_files, fasta)

    quant = lfq(
        mbr_fdr,
        fasta,
        as_sorted_files(mzml_files),
        as_sorted_files(id_files)
    )

    emit:
    id: Value<List<Path>> = as_sorted_files(mzml_files.mix(id_files))
    quant: Channel<Record> = quant
}
