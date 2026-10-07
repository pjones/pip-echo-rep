nextflow.enable.types = true

include { toMzML } from "./modules/thermo-to-mzml.nf"
include { flashlfq } from "./modules/flashlfq/main.nf"
include { proteomicslfq } from "./modules/proteomicslfq/main.nf"

params {
    dataset: String = "PXD037527"
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

    mzml_files = toMzML(mzml_ch)
    flfq = flashlfq(mzml_files.collect().map { fs -> fs.toSorted() }, fasta)
    plfq = proteomicslfq(mzml_files, fasta)

    publish:
    mzml_files = mzml_files
    flashlfq_before = flfq.quant_before
    flashlfq_after = flfq.quant_after
    proteomicslfq_before = plfq
}

output {
    mzml_files {
        path "${params.dataset}/mzML"
    }

    flashlfq_before {
        path "${params.dataset}/flashlfq/regular"
    }

    flashlfq_after {
        path "${params.dataset}/flashlfq/regular"
    }

    proteomicslfq_before {
        path "${params.dataset}/proteomicslfq/regular"
    }
}
