process PYSCENIC_CTX {
    tag "$sample"
    label 'process_high'
    container 'my-pyscenic:custom'
    publishDir "${params.outdir}/scenic/ctx", mode: 'copy'

    input:
    tuple val(sample), path(adj), path(h5ad), path(motif), path(annotation)

    output:
    tuple val(sample), path("${sample}_regulons.csv"), emit: regulons

    script:
    """
    pyscenic ctx ${adj} \\
        ${motif} \\
        --annotations_fname ${annotation} \\
        --expression_mtx_fname ${h5ad} \\
        --output ${sample}_regulons.csv \\
        --num_workers ${task.cpus}
    """
}