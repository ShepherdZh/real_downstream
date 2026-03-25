process CELLCHAT {
    tag "$sample_id"
    publishDir "${params.outdir}/cellchat", mode: 'copy'
    container "rocker/tidyverse:latest"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}.cellchat.rds")

    script:
    """
    run_cellchat.R \\
        --input ${h5ad} \\
        --output ${sample_id}.cellchat.rds \\
        --species human \\
        --cell_type_col cell_type
    """
}