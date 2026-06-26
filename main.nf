// Import your freshly tested subworkflow
include { VIDRL_UMI_DEDUP } from './subworkflows/local/vidrl_umi_dedup/main'

workflow {
    // 1. Create the paired FASTQ channel: [ val(meta), [ read_1, read_2 ] ]
    ch_fastq = Channel.fromFilePairs( params.input, checkIfExists: true )
        .map { id, reads ->
            def meta = [ id: id, single_end: false ]
            return [ meta, reads ]
        }

    // 2. Create the reference FASTA channel: [ val(meta_fasta), fasta_file ]
    ch_fasta = Channel.fromPath( params.fasta, checkIfExists: true )
        .map { fasta ->
            def meta_fasta = [ id: fasta.baseName ]
            return [ meta_fasta, fasta ]
        }

    // 3. Call your subworkflow with both channels
    VIDRL_UMI_DEDUP ( 
        ch_fastq, 
        ch_fasta 
    )
}