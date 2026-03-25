include { SCVELO as SCVELO_PROCESS } from '../modules/scvelo/main'

workflow SCVELO {
    take:
//    ch_input  // tuple val(sample), path(h5ad)

    main:
//    SCVELO_PROCESS(ch_input)
    SCVELO_PROCESS()
}
