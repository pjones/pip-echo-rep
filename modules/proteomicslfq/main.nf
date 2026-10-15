nextflow.enable.types = true

include { fasta_with_decoys } from "./decoy.nf"
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

workflow proteomicslfq {
    take:
    mbr_fdr: Float
    mzml_files: Channel<Path>
    fasta_files: List<Path>

    main:
    fasta = fasta_with_decoys(fasta_files)

    msgf_id = msgfplus(mzml_files, fasta)
    msgf_feat = msgfplus_features(msgf_id)

    comet_id = comet(mzml_files, fasta)
    comet_feat = comet_features(comet_id)

    def id_files_ch =
        msgf_feat
        .map { r -> record(base: r.base, msgfplus: r.percolator) }
        .join(comet_feat.map { r -> record(base: r.base, comet: r.percolator) }, by: "base")

    consensus_id = consensus(id_files_ch)
    id_before = mzml_files.mix(consensus_id)

    quant = lfq(
        mbr_fdr,
        as_sorted_files(mzml_files),
        as_sorted_files(consensus_id),
        fasta
    )

    emit:
    id_before: Channel<Path> = id_before
    quant_before: Channel<Record> = quant
}
