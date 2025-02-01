include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF   } from '../../modules/nf-core/gunzip'
include { CUSTOM_GTFFILTER       } from '../../modules/nf-core/custom/gtffilter'
include { STAR_GENOMEGENERATE    } from '../../modules/nf-core/star/genomegenerate'

workflow PREPARE_GENOME {
    main:
    ch_versions = Channel.empty()

    ch_fasta             = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    ch_gtf               = Channel.value([[id: 'gtf'], file(params.gtf, checkIfExists: true)])

    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( ch_fasta ).gunzip.collect()
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    }

    if (params.gtf.endsWith('.gz')) {
        ch_gtf      = GUNZIP_GTF ( ch_gtf ).gunzip.collect()
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    }

    CUSTOM_GTFFILTER (ch_gtf, ch_fasta)
    ch_gtf = CUSTOM_GTFFILTER.out.gtf.collect()
    ch_versions = ch_versions.mix(CUSTOM_GTFFILTER.out.versions)

    if (!params.star_index) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.star_index.collect()
    } else {
        ch_star_index = Channel.value([[id: 'star_index'], file(params.star_index, checkIfExists: true)])
    }

    emit:
    star_index = ch_star_index
    fasta      = ch_fasta
    gtf        = ch_gtf
    versions   = ch_versions
}
