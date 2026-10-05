include { toMzML } from "./modules/thermo-to-mzml.nf"

include {
    analysis as regular
} from "./modules/analysis.nf"

params {
    dataset: String = "PXD037527"
}

workflow {
    main:
    def spectra: Path = file("data/${params.dataset}/spectra.csv")

    def fasta: List<Path> = file("data/${params.dataset}/fasta.csv")
        .readLines()
        .collect { line ->
            file(spectra.parent.resolve(line))
        }

    def mzml_ch = channel.fromPath(spectra)
        .splitCsv()
        .map { line -> file(spectra.parent.resolve(line[0])) }

    toMzML(mzml_ch)
    regular(toMzML.out, fasta)

    publish:
    mzml_files = toMzML.out

    flashlfq_id_regular = regular.out.flashlfq_id
    flashlfq_quant_regular = regular.out.flashlfq_quant

    proteomicslfq_regular = regular.out.proteomicslfq
}

output {
    mzml_files {
        path "${params.dataset}/mzML"
    }

    flashlfq_id_regular {
        path "${params.dataset}/flashlfq/regular"
    }

    flashlfq_quant_regular {
        path "${params.dataset}/flashlfq/regular"
    }

    proteomicslfq_regular {
        path "${params.dataset}/proteomicslfq/regular"
    }
}
