process ANNDATA_CONCAT {
    tag "${meta.id}"

    label 'process_medium'

    conda "conda-forge::scanpy==1.10.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy:1.10.4--c2d474f46255931c' :
        'community.wave.seqera.io/library/scanpy:1.10.4--f905699eb17b6536' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'concat.py'

    stub:
    """
    touch combined_matrix.h5ad
    touch versions.yml
    """
}
