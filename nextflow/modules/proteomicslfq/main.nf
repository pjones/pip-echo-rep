include { fasta_with_decoys } from "./decoy.nf"
include { msgfplus } from "./msgfplus.nf"
include { comet } from "./comet.nf"
include { consensus } from "./consensus.nf"
include { proteomicslfq as lfq } from "./proteomicslfq.nf"

include {
    features as msgfplus_features ;
    features as comet_features
} from "./features.nf"

// zip :: Channel<Tuple<ID, a>, Tuple<ID, b>> -> Channel<Tuple<a, b>>
def zip(a, b) {
    a
        .map { t -> record(id: t[0], a: t[1]) }
        .join(b.map { t -> record(id: t[0], b: t[1]) }, by: 'id')
        .map { r -> [r.a, r.b] }
}

workflow proteomicslfq {
    take:
    mzml_files
    fasta_files

    main:
    fasta_with_decoys(fasta_files)

    msgfplus(mzml_files, fasta_with_decoys.out.database)
    msgfplus_features(msgfplus.out)

    comet(mzml_files, fasta_with_decoys.out.database)
    comet_features(comet.out)

    consensus(zip(msgfplus_features.out, comet_features.out))

    lfq(
        mzml_files.collect(),
        consensus.out.idxml.collect(),
        fasta_with_decoys.out.database,
    )

    emit:
    lfq.out
}
