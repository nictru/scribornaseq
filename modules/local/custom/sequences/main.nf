process CUSTOM_SEQUENCES {
    tag "${meta.id}"
    label 'process_low'

    conda 'conda-forge::biopython=1.85 conda-forge::pyyaml=6.0.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_pyyaml:9805990d7ff97d63' :
        'community.wave.seqera.io/library/biopython_pyyaml:1e6c45b1935e2dd4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.gtf")  , emit: gtf
    path "versions.yml"                     , emit: versions

    script:
    prefix = task.ext.prefix ?: "$meta.id"
    template "sequences.py"
}
