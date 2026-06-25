// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FASTP     } from '../../modules/nf-core/fastp/main'
include { BWA_MEM   } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { FGBIO_EXTRACTUMISFROMREADSTRUCTURE } from '../../modules/nf-core/fgbio/extractumisfromreadstructure/main'
include { FGBIO_GROUPREADSBYUMI            } from '../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main'

workflow VIDRL_UMI_DEDUP {

    take:
    ch_fastq // channel: [ val(meta), [ read_1, read_2 ] ]
    ch_fasta // channel: [ val(meta_fasta), fasta_file ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow
    ch_versions = Channel.empty()

    // Step 1: Generate BWA index on-the-fly from the reference FASTA
    BWA_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    // Step 2: Quality trimming on raw paired-end reads
    FASTP ( ch_fastq, [], false, false )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // Step 3: Combine trimmed reads with the dynamically generated index
    // This cross-references every sample meta with the same index file
    ch_bwa_in = FASTP.out.reads.combine(BWA_INDEX.out.index)

    // Step 4: Align reads using the new index 
    BWA_MEM (
        ch_bwa_in,
        ch_fasta,
        true // sort_bam: outputs coordinate-sorted BAM needed for fgbio
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // 5. Extract UMI from Read ID and write into the SAM 'RX' Tag
    // fgbio requires a read-structure config (e.g., via modules.config)
    FGBIO_EXTRACTUMISFROMREADSTRUCTURE ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(FGBIO_EXTRACTUMISFROMREADSTRUCTURE.out.versions)

    // 6. Group aligned, tagged reads by coordinate + UMI similarity
    FGBIO_GROUPREADSBYUMI ( FGBIO_EXTRACTUMISFROMREADSTRUCTURE.out.bam )
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions)

    //7. Strip duplicates completely using UMI information embedded in the groups
    PICARD_MARKDUPLICATES ( FGBIO_GROUPREADSBYUMI.out.bam, [[:],[]], [[:],[]] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    bam      = PICARD_MARKDUPLICATES.out.bam      // Grouped BAM ready for deduplication
    metrics  = FGBIO_GROUPREADSBYUMI.out.histogram // UMI family size distribution metrics
    versions = ch_versions                        // Software versions channel
}
