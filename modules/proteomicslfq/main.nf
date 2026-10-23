nextflow.enable.types = true

include { fasta_with_decoys } from "./decoy.nf"
include { censor } from "./censor.nf"
include { analysis } from "../analysis.nf"

include {
    id_quant as id_quant_before ;
    id_quant as id_quant_after
} from "./id_quant.nf"

workflow proteomicslfq {
    take:
    mbr_fdr: Channel<Float>
    mzml_files: Channel<Path>
    fasta_files: List<Path>

    main:
    fasta = fasta_with_decoys(fasta_files)
    before = id_quant_before(mbr_fdr, fasta, mzml_files)

    censored = censor(before.id)
    censored_tsv = censored.map {r -> r.censored}
    censored_mzml = censored.flatMap {r -> r.mzml}

    after = id_quant_after(mbr_fdr, fasta, censored_mzml)

    before_quant = before.quant
        .map { r ->
            record(
                mbr_fdr: r.mbr_fdr,
                orig_path: r.consensus.parent
            )
        }

    after_quant = after.quant.map {r ->
        record(
            mbr_fdr: r.mbr_fdr,
            censored_path: r.consensus.parent
        )
    }

    analyzed = analysis(
        "proteomicslfq",
        params.data_dir,
        censored.map { r -> r.censored },
        before_quant.join(after_quant, by: "mbr_fdr")
    )

    emit:
    before_id = before.id
    before_quant = before.quant
    censored = censored_tsv
    after_id = after.id
    after_quant = after.quant
    analyzed: Channel<Record> = analyzed
}
