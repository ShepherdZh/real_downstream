process SCVELO {
    tag "null"
    container "my-scvelo:latest"
    publishDir "${params.outdir}/scvelo", mode: 'copy'

    script:
    """
    run_scvelo.py \\
        --output sample_velocity.h5ad \\
        --mode dynamical
    """
}