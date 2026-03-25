include { PYSCENIC_GRN }   from '../modules/scenic/pyscenic_grn'
include { PYSCENIC_CTX }   from '../modules/scenic/pyscenic_ctx'
include { PYSCENIC_AUCELL } from '../modules/scenic/pyscenic_aucell'

workflow SCENIC {
    take:
    h5ad_files  // tuple val(sample), path(h5ad)
    tfs_ch    // path to TF list file
    motif_db_ch
    annotation_db_ch
	
    main:
    // Step 1: GRN inference - combine h5ad with TF file
    PYSCENIC_GRN(h5ad_files.combine(tfs_ch))
    
    // Step 2: Motif enrichment (cisTarget)
    def ch_ctx_input = PYSCENIC_GRN.out.adjacencies
        .join(h5ad_files, by: 0)  // 按 sample (索引0) join
        .combine(motif_db_ch)
        .combine(annotation_db_ch)
    // 现在顺序是: [sample, adj.tsv, h5ad, motif_db, annotation_db]。貌似adj和h5ad都是带着sample id的tuple，这俩不能用combine，否则sample id出现两次；join 将基于 sample id 合并两个通道
    
    PYSCENIC_CTX(ch_ctx_input)
    
    // Step 3: AUCell
    def ch_for_auc = h5ad_files.join(PYSCENIC_CTX.out.regulons, by: 0)
    PYSCENIC_AUCELL(ch_for_auc)

//不做emit:，全交给module的process去emit了
}