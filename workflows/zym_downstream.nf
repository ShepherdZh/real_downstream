//include { CCC      } from '../subworkflows/ccc'
include { SCENIC   } from '../subworkflows/scenic'
include { SCVELO   } from '../subworkflows/scvelo'

workflow SCRNA_DOWNSTREAM {

    take:
    h5ad_files

    main:
    // Create TF channel from params
    def tfs_ch = channel.fromPath(params.scenic_tfs, checkIfExists: true)
    def annotation_db_ch = channel.fromPath(params.scenic_annotation_db, checkIfExists: true)
    def motif_db_ch = channel.fromPath(params.scenic_motif_db, checkIfExists: true)
    
//    CCC(h5ad_files)
    SCENIC(h5ad_files, tfs_ch,motif_db_ch,annotation_db_ch)
//    SCVELO(h5ad_files)
    SCVELO()
}