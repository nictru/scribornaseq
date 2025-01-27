process BAM_PRIORITIZE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_polars_pyyaml:b25c7c0c4fc463b6' :
        'community.wave.seqera.io/library/pysam_polars_pyyaml:a89a5cd920894783' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)
    val gene_type

    output:
    tuple val(meta), path("${prefix}.modified.bam") , emit: modified
    tuple val(meta), path("${prefix}.colliders.bam"), emit: colliders
    path "versions.yml"                             , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "prioritize.py"
}
