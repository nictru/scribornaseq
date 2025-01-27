process BAM_MULTIMAPPERS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hdcf5f25_2' :
        'biocontainers/bamtools:2.5.2--hdcf5f25_2' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.multimappers.bam"), emit: multimappers
    tuple val(meta), path("${prefix}.unique.bam"),       emit: uniquemappers
    path  "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamtools filter -in ${bam} -mapQuality '=255' -out ${prefix}.unique.bam
    bamtools filter -in ${bam} -mapQuality '<255' -out ${prefix}.multimappers.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
