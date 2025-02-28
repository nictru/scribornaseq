include { STAR_STARSOLO as FIRST_PASS  } from '../../modules/local/star/starsolo'
include { STAR_STARSOLO as SECOND_PASS } from '../../modules/local/star/starsolo'

workflow STAR2PASS {

    take:
    ch_reads
    ch_star_index
    ch_gtf
    barcode_whitelist
    protocol
    star_feature

    main:
    ch_versions = Channel.empty()

    FIRST_PASS(
        ch_reads,
        ch_star_index,
        ch_gtf,
        barcode_whitelist,
        protocol,
        star_feature,
        []
    )
    ch_versions = ch_versions.mix(FIRST_PASS.out.versions)

    SECOND_PASS(
        ch_reads,
        ch_star_index,
        ch_gtf,
        barcode_whitelist,
        protocol,
        star_feature,
        FIRST_PASS.out.sjdb.map{it[1]}.collect()
    )
    ch_versions = ch_versions.mix(SECOND_PASS.out.versions)

    emit:
    raw_counts = SECOND_PASS.out.raw_counts
    raw_velocyto = SECOND_PASS.out.raw_velocyto
    raw_sj = SECOND_PASS.out.raw_sj
    versions = ch_versions
}