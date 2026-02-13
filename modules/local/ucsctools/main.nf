

process FA_TO_TWOBIT {
    tag "$fasta"
    label 'process_low'

    /*
    container not available?
    */
    conda "${moduleDir}/environment.yml"

    input:
    path fasta

    output:
    path "*.2bit", emit: twobit
    path "versions.yml", emit: versions

    script:
    def prefix = fasta.baseName
    """
    faToTwoBit $fasta ${prefix}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faToTwoBit: \$(faToTwoBit 2>&1 | head -n 1 | sed 's/faToTwoBit - //')
    END_VERSIONS
    """

    stub:
    def prefix = fasta.baseName
    """
    touch ${prefix}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        faToTwoBit: \$(faToTwoBit 2>&1 | head -n 1 | sed 's/faToTwoBit - //')
    END_VERSIONS
    """
}



