process PYSCENIC_GRN {
    tag "$sample"
    label 'process_high'
    container 'my-pyscenic:custom'
    publishDir "${params.outdir}/scenic/grn", mode: 'copy'

    input:
    tuple val(sample), path(h5ad), path(tfs)//次序要与之前include的函数的参数顺序一致

    output:
    tuple val(sample), path("${sample}_adj.tsv"), emit: adjacencies

    script:
    """
    pyscenic grn ${h5ad} \\
        ${tfs} \\
        --num_workers ${task.cpus} \\
        --output ${sample}_adj.tsv \\
        --seed 123
    """
}