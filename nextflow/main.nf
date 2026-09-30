include { toMzML } from "./modules/thermo-to-mzml.nf"
include { flashlfq } from "./modules/flashlfq/main.nf"
include { proteomicslfq } from "./modules/proteomicslfq/main.nf"

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

    mzml_ch = channel.fromPath(spectra)
        .splitCsv()
        .map { line -> file(spectra.parent.resolve(line[0])) }

    toMzML(mzml_ch)
    flashlfq(toMzML.out, fasta)
    plfq_res = proteomicslfq(toMzML.out, fasta)

    publish:
    mzml_files = toMzML.out
    mztab = plfq_res.mztab
}

output {
    mzml_files {
        path "${params.dataset}/mzML"
    }
    mztab {
        path "${params.dataset}/proteomicslfq.mzTab"
    }
}
