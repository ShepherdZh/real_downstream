include { CELLCHAT } from '../modules/cellchat/main'

workflow CCC {
    take:
    ch_input  // tuple val(sample), path(h5ad)

    main:
    CELLCHAT(ch_input)

    emit:
    cellchat_rds = CELLCHAT.out
}
