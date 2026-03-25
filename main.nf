include { SCRNA_DOWNSTREAM } from './workflows/zym_downstream'

workflow {
    // Create input channel from h5ad files
    def h5ad_files = channel
        .fromPath(params.input, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    // Run downstream analysis workflow
    SCRNA_DOWNSTREAM(h5ad_files)
}
