include { FASTQC } from '../../modules/nf-core/fastqc'
include { FASTP  } from '../../modules/nf-core/fastp'


workflow PREPARE_READS {
    take:
    ch_reads

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (!params.skip_fastqc) {
        FASTQC (
            ch_reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if (!params.skip_fastp) {
        ch_split_reads = ch_reads.multiMap{meta, reads ->
            read1: [meta, reads[0]]
            read2: [meta + [single_end: true], reads[1]]
        }
        FASTP (
            ch_split_reads.read2,
            [], [], [], []
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
        ch_versions = ch_versions.mix(FASTP.out.versions)

        ch_reads = ch_split_reads.read1.join(FASTP.out.reads.map{meta, read -> [meta + [single_end: false], read]})
            .map{meta, read1, read2 -> [meta, [read1, read2]]}
    }

    emit:
    reads = ch_reads
    versions = ch_versions
    multiqc_files = ch_multiqc_files
}
