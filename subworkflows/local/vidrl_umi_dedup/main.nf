// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
include { SAMTOOLS_FAIDX     } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { FASTP     } from '../../../modules/nf-core/fastp/main'
include { BWA_MEM   } from '../../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main'
include { FGBIO_COPYUMIFROMREADNAME } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_PRE_GROUP  } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_POST_GROUP } from '../../../modules/nf-core/samtools/sort/main'
include { FGBIO_GROUPREADSBYUMI            } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'

workflow VIDRL_UMI_DEDUP {

    take:
    ch_fastq // channel: [ val(meta), [ read_1, read_2 ] ]
    ch_fasta // channel: [ val(meta_fasta), fasta_file ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow
    //ch_versions = Channel.empty()

    // Step 1: Generate BWA index on-the-fly from the reference FASTA
    BWA_INDEX ( ch_fasta )
    //ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    ch_fasta_raw = ch_fasta.map{
        meta, fasta -> fasta
    }

    ch_fai = SAMTOOLS_FAIDX(
        ch_fasta.map { meta, fasta -> [ [id: fasta.baseName], fasta, [] ] },  
        false
    ).fai

    // Step 2: Quality trimming on raw paired-end reads
    //FASTP ( ch_fastq, [], false, false )
    FASTP(
        ch_fastq.map{ meta, reads ->
            [meta, reads, []]
        },
        false,
        false,
        false
    )
    //ch_versions = ch_versions.mix(FASTP.out.versions)

    // Step 3: Combine trimmed reads with the dynamically generated index
    // This cross-references every sample meta with the same index file
    //ch_bwa_in = FASTP.out.reads.combine(BWA_INDEX.out.index)

    // Step 4: Align reads using the new index 
    BWA_MEM (
        FASTP.out.reads,
        BWA_INDEX.out.index,
        ch_fasta,
        true // sort_bam: outputs coordinate-sorted BAM needed for fgbio
    )
    SAMTOOLS_INDEX( BWA_MEM.out.bam )
    //ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    ch_bam_with_index = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.index)
    // 5. Extract UMI from Read ID and write into the SAM 'RX' Tag
    // fgbio requires a read-structure config (e.g., via modules.config)
    FGBIO_COPYUMIFROMREADNAME ( ch_bam_with_index )
    //ch_versions = ch_versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)

    SAMTOOLS_SORT_PRE_GROUP ( 
        FGBIO_COPYUMIFROMREADNAME.out.bam, 
        [ [:], [], [] ],               // Matches input 2: Empty reference tuple
        'bai' 
    )
    // 6. Group aligned, tagged reads by coordinate + UMI similarity
    FGBIO_GROUPREADSBYUMI ( 
        SAMTOOLS_SORT_PRE_GROUP.out.bam,
        "paired" //paired strategy is selected
    )

    SAMTOOLS_SORT_POST_GROUP ( 
        FGBIO_GROUPREADSBYUMI.out.bam, 
        [ [:], [], [] ],               // Matches input 2: Empty reference tuple
        'bai' 
    )
    //ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions)

    picard_reference = ch_fasta
    .combine(ch_fai)
    //.map { it ->                    // it is the full tuple [fasta, fai]
    //    def fasta = it[0]
    //    def fai   = it[1]
    //   [ [id: fasta.baseName], fasta, fai ]
    .map { meta_fasta, fasta, meta_fai, fai ->
        // Explicitly structural alignment: [ meta, fasta, fai ]
        return [ meta_fasta, fasta, fai ]
    }
    .collect()
    //7. Strip duplicates completely using UMI information embedded in the groups
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT_POST_GROUP.out.bam, picard_reference)
    //ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    bam      = PICARD_MARKDUPLICATES.out.bam      // Grouped BAM ready for deduplication
    metrics  = FGBIO_GROUPREADSBYUMI.out.histogram // UMI family size distribution metrics
    //versions = ch_versions                        // Software versions channel
}
