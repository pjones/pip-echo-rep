nextflow.enable.types = true

process flashlfq {
    input:
    mbr_fdr: Float
    record(
      t3_psms: Path,
      t2_mzml: Set<Path>,
    )

    output:
    record(
        mbr_fdr: mbr_fdr,
        peaks: file("QuantifiedPeaks.tsv"),
        peptides: file("QuantifiedPeptides.tsv"),
        proteins: file("QuantifiedProteins.tsv")
    )

    script:
    // NOTE: The mzML files task 2 (t2_mzml) will appear in the
    // current directory as symbolic links.  Therefore they are not
    // directly referenced below but instead the current directory is
    // given to the `-d` option.
    """
    flashlfq.sh \
        -d . \
        -p "${t3_psms}" \
        -F "${mbr_fdr}" \
        -I "${params.isotope_ppm_tolerance}" \
        -N "${params.num_isotopes_required}" \
        -P "${params.ppm_tolerance}" \
        -W "${params.mbr_rt_window}"
    """
}
