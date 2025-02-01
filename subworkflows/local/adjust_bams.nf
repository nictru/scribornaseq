include { SAMTOOLS_VIEW as EXTRACT_MULTIMAPPERS  } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW as EXTRACT_UNIQUEMAPPERS } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX                         } from '../../modules/nf-core/samtools/index'
include { BAM_PRIORITIZE                         } from '../../modules/local/bam/prioritize'
include { SAMTOOLS_MERGE                         } from '../../modules/nf-core/samtools/merge'
include { SAMTOOLS_SORT                          } from '../../modules/nf-core/samtools/sort'
include { PICARD_CLEANSAM                        } from '../../modules/nf-core/picard/cleansam'


workflow ADJUST_BAMS {
    take:
    ch_bam
    ch_gtf
    ch_fasta

    main:
    ch_versions = Channel.empty()

    ch_bam_bai = ch_bam.map{meta, bam -> [meta, bam, []]}

    EXTRACT_MULTIMAPPERS(ch_bam_bai, [[], []], [])
    ch_versions = ch_versions.mix(EXTRACT_MULTIMAPPERS.out.versions)

    EXTRACT_UNIQUEMAPPERS(ch_bam_bai, [[], []], [])
    ch_versions = ch_versions.mix(EXTRACT_UNIQUEMAPPERS.out.versions)

    ch_multimappers  = EXTRACT_MULTIMAPPERS.out.bam
    ch_uniquemappers = EXTRACT_UNIQUEMAPPERS.out.bam

    BAM_PRIORITIZE(ch_multimappers.join(EXTRACT_MULTIMAPPERS.out.csi), ch_gtf, "rRNA")
    ch_versions = ch_versions.mix(BAM_PRIORITIZE.out.versions)

    SAMTOOLS_SORT(BAM_PRIORITIZE.out.modified, ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    ch_bams = ch_uniquemappers.join(SAMTOOLS_SORT.out.bam)
        .map{meta, unique, multi -> [meta, [unique, multi]]}

    SAMTOOLS_MERGE(ch_bams, [[], []], [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    PICARD_CLEANSAM(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(PICARD_CLEANSAM.out.versions)


    emit:
    bam      = PICARD_CLEANSAM.out.bam
    versions = ch_versions
}
