process LOESS_SCORE {
    tag "$meta_sample.sampleid"
    label 'process_medium'

    input:
    tuple val(meta_sample), val(source), path(matrix)

    output:
    tuple val(meta_sample), val(source), path("*_loess_score.json"), emit: loess_score_json
    tuple val(meta_sample), val(source), path("*_spans.pdf"), path("*_loess.pdf"), emit: loess_score_pdfs
    tuple val(meta_sample), val(source), path("*_spans.png"), path("*_loess.png"), emit: loess_score_pngs
    path "versions.yml", emit: versions

    script:
    def span_arg = task.ext.span ? "--span ${task.ext.span}" : ''
    """
    for MAT in ${matrix.join(' ')}; do
        Rscript ${moduleDir}/resources/usr/bin/loess_score.r \\
            --file \$MAT \\
            ${span_arg}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version | sed -e "s/Rscript (R) //g")
    END_VERSIONS
    """

    stub:
    """
    for MAT in ${matrix.join(' ')}; do
        BASENAME=\$(basename \$MAT _matrix.gz)
        OUTPUT_STUB=\${BASENAME}_\${meta_sample.sampleid}
        touch \${OUTPUT_STUB}_loess_score.json
        touch \${OUTPUT_STUB}_spans.pdf
        touch \${OUTPUT_STUB}_loess.pdf
        touch \${OUTPUT_STUB}_spans.png
        touch \${OUTPUT_STUB}_loess.png
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version | sed -e "s/Rscript (R) //g")
    END_VERSIONS
    """
}
