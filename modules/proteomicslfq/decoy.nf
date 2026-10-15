nextflow.enable.types = true

// Create a protein database with decoys.
process fasta_with_decoys {
    input:
    fasta_files: List<Path>

    output:
    file("database.revCat.fasta")

    script:
    """
    msgfplus=\$(which msgfplus)
    jar=\$(realpath "\$(dirname "\$msgfplus")/../share/msgfplus/MSGFPlus.jar")

    cat ${fasta_files.join(' ')} >"database.fasta"

    java \
      -Xms6144m \
      -Xmx6144m \
      -cp "\$jar" \
      edu.ucsd.msjava.msdbsearch.BuildSA \
      -d "database.fasta" \
      -o . \
      -debug 0 \
      -tda 1 \
      -decoy "DECOY_"
    """
}
