nextflow.enable.types = true

include { fasta_with_decoys } from "./decoy.nf"
include { censor } from "./censor.nf"

include {
    id_quant as id_quant_before ;
    id_quant as id_quant_after
} from "./id_quant.nf"

workflow proteomicslfq {
    take:
    mbr_fdr: Float
    mzml_files: Channel<Path>
    fasta_files: List<Path>

    main:
    fasta = fasta_with_decoys(fasta_files)
    before = id_quant_before(mbr_fdr, fasta, mzml_files)

    censored = censor(before.id)
    censored_mzml = censored.flatMap {r -> r.mzml}

    after = id_quant_after(mbr_fdr, fasta, censored_mzml)

    emit:
    before = before
    censored = censored
    after = after
}
