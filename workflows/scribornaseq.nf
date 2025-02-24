/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_READS                         } from '../subworkflows/local/prepare_reads'
include { PREPARE_GENOME                        } from '../subworkflows/local/prepare_genome'
include { STAR_STARSOLO                         } from '../modules/local/star/starsolo'
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

    ch_barcode_whitelist = Channel.value(file(params.barcode_whitelist, checkIfExists: true))

    PREPARE_READS(ch_samplesheet)
    ch_reads = PREPARE_READS.out.reads
    ch_versions = ch_versions.mix(PREPARE_READS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_READS.out.multiqc_files)

    PREPARE_GENOME()
    ch_star_index = PREPARE_GENOME.out.star_index
    ch_fasta      = PREPARE_GENOME.out.fasta
    ch_gtf        = PREPARE_GENOME.out.gtf
    ch_versions   = ch_versions.mix(PREPARE_GENOME.out.versions)

    STAR_STARSOLO(
        ch_reads,
        ch_star_index,
        ch_gtf,
        ch_barcode_whitelist,
        "CB_UMI_Simple",
        params.solo_feature
    )
    ch_versions = ch_versions.mix(STAR_STARSOLO.out.versions)

    ch_counts = STAR_STARSOLO.out.raw_counts
        .join(STAR_STARSOLO.out.raw_velocyto, remainder: true)
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
