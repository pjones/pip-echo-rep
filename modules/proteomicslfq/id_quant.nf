nextflow.enable.types = true

include { msgfplus } from "./msgfplus.nf"
include { comet } from "./comet.nf"
include { consensus } from "./consensus.nf"
include { proteomicslfq as lfq } from "./proteomicslfq.nf"

include {
    features as msgfplus_features ;
    features as comet_features
} from "./features.nf"

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
    msgf_id = msgfplus(mzml_files, fasta)
    msgf_feat = msgfplus_features(msgf_id)

    comet_id = comet(mzml_files, fasta)
    comet_feat = comet_features(comet_id)

    def id_files_ch =
        msgf_feat
        .map { r -> record(base: r.base, msgfplus: r.percolator) }
        .join(comet_feat.map { r -> record(base: r.base, comet: r.percolator) }, by: "base")

    consensus_id = consensus(id_files_ch)
    id = as_sorted_files(mzml_files.mix(consensus_id))

    quant = lfq(
        mbr_fdr,
        as_sorted_files(mzml_files),
        as_sorted_files(consensus_id),
        fasta
    )

    emit:
    id: Value<List<Path>> = id
    quant: Channel<Record> = quant
}
