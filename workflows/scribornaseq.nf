/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                } from '../modules/nf-core/fastqc'
include { FASTP                                 } from '../modules/nf-core/fastp'
include { GUNZIP as GUNZIP_FASTA                } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                  } from '../modules/nf-core/gunzip'
include { CUSTOM_GTFFILTER                      } from '../modules/nf-core/custom/gtffilter'
include { STAR_GENOMEGENERATE                   } from '../modules/nf-core/star/genomegenerate'
include { STAR_STARSOLO as ALIGN                } from '../modules/local/star/starsolo'
include { SAMTOOLS_VIEW as EXTRACT_MULTIMAPPERS } from '../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW as EXTRACT_UNIQUEMAPPERS} from '../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX                        } from '../modules/nf-core/samtools/index'
include { BAM_PRIORITIZE                        } from '../modules/local/bam/prioritize'
include { SAMTOOLS_MERGE                        } from '../modules/nf-core/samtools/merge'
include { SAMTOOLS_SORT                         } from '../modules/nf-core/samtools/sort'
include { PICARD_CLEANSAM                       } from '../modules/nf-core/picard/cleansam'
include { STAR_STARSOLO as QUANTIFY             } from '../modules/local/star/starsolo'
include { ANNDATA_READMTX                       } from '../modules/local/anndata/readmtx'
include { ANNDATA_CONCAT                        } from '../modules/local/anndata/concat'
include { MULTIQC                               } from '../modules/nf-core/multiqc'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_scribornaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCRIBORNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    ch_reads = ch_samplesheet.multiMap{meta, reads ->
        read1: [meta, reads[0]]
        read2: [meta + [single_end: true], reads[1]]
    }
    FASTP (
        ch_reads.read2,
        [], [], [], []
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    ch_trimmed = ch_reads.read1.join(FASTP.out.reads.map{meta, read -> [meta + [single_end: false], read]})
        .map{meta, read1, read2 -> [meta, [read1, read2]]}

    ch_fasta             = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)])
    ch_gtf               = Channel.value([[id: 'gtf'], file(params.gtf, checkIfExists: true)])
    ch_barcode_whitelist = Channel.value(file(params.barcode_whitelist, checkIfExists: true))

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

    ALIGN(
        ch_trimmed,
        ch_star_index,
        ch_gtf,
        ch_barcode_whitelist,
        "CB_UMI_Simple",
        params.solo_feature
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    // Samttols view expects an index
    ch_alignment_bam = ALIGN.out.bam_sorted.map{meta, bam -> [meta, bam, []]}

    EXTRACT_MULTIMAPPERS(ch_alignment_bam, [[], []], [])
    ch_versions = ch_versions.mix(EXTRACT_MULTIMAPPERS.out.versions)

    EXTRACT_UNIQUEMAPPERS(ch_alignment_bam, [[], []], [])
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

    QUANTIFY(
        PICARD_CLEANSAM.out.bam,
        ch_star_index,
        ch_gtf,
        ch_barcode_whitelist,
        "CB_UMI_Simple",
        params.solo_feature
    )
    ch_versions = ch_versions.mix(QUANTIFY.out.versions)

    ch_counts = QUANTIFY.out.raw_counts
        .join(QUANTIFY.out.raw_velocyto, remainder: true)
        .map{meta, count, velocity -> {
                return [meta + [input_type: 'raw'], velocity ? [count, velocity] : [count]]
            }
        }

    ANNDATA_READMTX(ch_counts)
    ch_versions = ch_versions.mix(ANNDATA_READMTX.out.versions)

    ANNDATA_CONCAT(ANNDATA_READMTX.out.h5ad.map{it[1]}.collect().map{[[id: 'combined'], it]})
    ch_versions = ch_versions.mix(ANNDATA_CONCAT.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'scribornaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
