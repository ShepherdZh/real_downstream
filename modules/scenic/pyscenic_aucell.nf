process PYSCENIC_AUCELL {
    tag "$sample"
    label 'process_medium'
    container 'my-pyscenic:custom'
    publishDir "${params.outdir}/scenic/aucell", mode: 'copy'

    input:
    tuple val(sample), path(h5ad), path(motifs)  // join 后顺序: sample, h5ad, motifs

    output:
    tuple val(sample), path("${sample}_auc_mtx.csv"), emit: auc

    script:
    """
    pyscenic aucell ${h5ad} \\
        ${motifs} \\
        --num_workers ${task.cpus} \\
        --output ${sample}_auc_mtx.csv
    """
}