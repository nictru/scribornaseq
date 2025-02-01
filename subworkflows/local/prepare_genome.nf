include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF          } from '../../modules/nf-core/gunzip'
include { CUSTOM_SEQUENCES              } from '../../modules/local/custom/sequences'
include { CAT_CAT as CONCAT_FASTA       } from '../../modules/nf-core/cat/cat'
include { CAT_CAT as CONCAT_SEQ_GTF     } from '../../modules/nf-core/cat/cat'
include { CAT_CAT as CONCAT_CUSTOM_GTF  } from '../../modules/nf-core/cat/cat'
include { GNU_SORT as SORT_GTF          } from '../../modules/nf-core/gnu/sort'
include { CUSTOM_GTFFILTER              } from '../../modules/nf-core/custom/gtffilter'
include { STAR_GENOMEGENERATE           } from '../../modules/nf-core/star/genomegenerate'

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

    if (params.custom_annotations) {
        ch_custom_gtf = Channel.value(file(params.custom_annotations, checkIfExists: true))
        CONCAT_CUSTOM_GTF(ch_gtf.combine(ch_custom_gtf)
            .map{meta, gtf, custom_gtf -> [meta, [gtf, custom_gtf]]})
        ch_gtf = CONCAT_CUSTOM_GTF.out.file_out.collect()
        ch_versions = ch_versions.mix(CONCAT_CUSTOM_GTF.out.versions)
    }

    if (params.custom_sequences) {
        CUSTOM_SEQUENCES (Channel.value([[id: 'custom_sequences'], file(params.custom_sequences, checkIfExists: true)]))
        ch_versions = ch_versions.mix(CUSTOM_SEQUENCES.out.versions)

        CONCAT_FASTA(ch_fasta.combine(CUSTOM_SEQUENCES.out.fasta.map{it[1]})
            .map{meta, genome_fasta, custom_fasta -> [meta, [genome_fasta, custom_fasta]]})
        ch_fasta = CONCAT_FASTA.out.file_out.collect()
        ch_versions = ch_versions.mix(CONCAT_FASTA.out.versions)

        CONCAT_SEQ_GTF(ch_gtf.combine(CUSTOM_SEQUENCES.out.gtf.map{it[1]})
            .map{meta, genome_gtf, custom_gtf -> [meta, [genome_gtf, custom_gtf]]})
        ch_gtf = CONCAT_SEQ_GTF.out.file_out.collect()
        ch_versions = ch_versions.mix(CONCAT_SEQ_GTF.out.versions)
    }

    if (params.custom_annotations || params.custom_sequences) {
        SORT_GTF(ch_gtf)
        ch_gtf = SORT_GTF.out.sorted.collect()
        ch_versions = ch_versions.mix(SORT_GTF.out.versions)
    }

    CUSTOM_GTFFILTER (ch_gtf, ch_fasta)
    ch_gtf = CUSTOM_GTFFILTER.out.gtf.collect()
    ch_versions = ch_versions.mix(CUSTOM_GTFFILTER.out.versions)

    if (!params.star_index) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_star_index = STAR_GENOMEGENERATE.out.index.collect()
    } else {
        ch_star_index = Channel.value([[id: 'star_index'], file(params.star_index, checkIfExists: true)])
    }

    emit:
    star_index = ch_star_index
    fasta      = ch_fasta
    gtf        = ch_gtf
    versions   = ch_versions
}
