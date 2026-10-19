nextflow.enable.types = true

include { toMzML } from "./modules/thermo-to-mzml.nf"
include { flashlfq } from "./modules/flashlfq/main.nf"
include { proteomicslfq } from "./modules/proteomicslfq/main.nf"

params {
    fdrs: List<Float> = [
        0.01,
        0.025,
        0.05,
        1.0, // MBR FDR disabled
    ]

    dataset: String = "PXD037527"
    data_dir: Path = file("data")
    human_file_pattern: String = "_1x02nguL_"
    foreign_species: String = "Saccharomyces cerevisiae"
    rt_delta: Float = 40 * 0.01

    ppm_tolerance: Float = 10.0
    isotope_ppm_tolerance: Float = 5.0
    num_isotopes_required: Integer = 2
    mbr_rt_window: Float = 1.5 // Minutes
}

workflow {
    main:
    def spectra: Path = file("data/${params.dataset}/spectra.csv")
    def mzml_dir: Path = spectra.parent

    def fasta = file("data/${params.dataset}/fasta.csv")
        .readLines()
        .collect { line -> mzml_dir.resolve(line) }

    def mzml_ch = channel.fromPath("${spectra}")
        .flatMap { csv -> csv.splitCsv(header: ["file"], skip: 0) }
        .map { row -> mzml_dir.resolve(row.file) }

    def fdr_ch = channel.fromList(params.fdrs)

    mzml_files = toMzML(mzml_ch)

    // MBR FDR cannot be disabled in FlashLFQ:
    def flash_fdr = fdr_ch.filter { fdr -> fdr < 1.0 }

    flfq = flashlfq(flash_fdr, mzml_files.collect().map { fs -> fs.toSorted() }, fasta)
    plfq = proteomicslfq(fdr_ch, mzml_files, fasta)

    publish:
    mzml_files = mzml_files

    flashlfq_analyzed = flfq.analyzed
    flashlfq_before_id = flfq.id_before
    flashlfq_before_quant = flfq.quant_before
    flashlfq_after_id = flfq.id_after
    flashlfq_after_quant = flfq.quant_after

    proteomicslfq_before_id = plfq.before_id
    proteomicslfq_before_quant = plfq.before_quant
    proteomicslfq_censored = plfq.censored
    proteomicslfq_after_id = plfq.after_id
    proteomicslfq_after_quant = plfq.after_quant
}

output {
    mzml_files {
        path "${params.dataset}/mzML"
    }

    flashlfq_before_id {
        path "${params.dataset}/flashlfq/original"
    }

    flashlfq_before_quant {
        path { quant -> "${params.dataset}/flashlfq/original/fdr-${quant.mbr_fdr}" }
    }

    flashlfq_after_id {
        path "${params.dataset}/flashlfq/censored"
    }

    flashlfq_after_quant {
        path { quant -> "${params.dataset}/flashlfq/censored/fdr-${quant.mbr_fdr}" }
    }

    flashlfq_analyzed {
        path "${params.dataset}/flashlfq/analyzed"
    }

    proteomicslfq_before_id {
        path "${params.dataset}/proteomicslfq/original/id"
    }

    proteomicslfq_before_quant {
        path { quant -> "${params.dataset}/proteomicslfq/original/fdr-${quant.mbr_fdr}" }
    }

    proteomicslfq_censored {
        path "${params.dataset}/proteomicslfq/censored/id"
    }

    proteomicslfq_after_id {
        path "${params.dataset}/proteomicslfq/censored/id"
    }

    proteomicslfq_after_quant {
        path { quant -> "${params.dataset}/proteomicslfq/censored/fdr-${quant.mbr_fdr}" }
    }
}
