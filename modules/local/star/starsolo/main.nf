process STAR_STARSOLO {

    //
    // This module executes STAR align quantification
    //

    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::star=2.7.11b bioconda::samtools=1.21'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_star:49284a612da8d71d' :
        'community.wave.seqera.io/library/samtools_star:bc8099e10cdd1c08' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    path whitelist
    val protocol
    val star_feature
    path sjdb

    output:
    tuple val(meta), path('*d.out.bam')                            , emit: bam
    tuple val(meta), path('*.Solo.out')                            , emit: counts
    tuple val(meta), path ("*.Solo.out/Gene*/raw")                 , emit: raw_counts
    tuple val(meta), path ("*.Solo.out/Gene*/filtered")            , emit: filtered_counts
    tuple val(meta), path ("*.Solo.out/SJ/sj_raw")                 , emit: raw_sj, optional:true
    tuple val(meta), path ("*.Solo.out/Velocyto/velocyto_raw")     , emit: raw_velocyto, optional:true
    tuple val(meta), path ("*.Solo.out/Velocyto/velocyto_filtered"), emit: filtered_velocyto, optional:true
    tuple val(meta), path('*Log.final.out')                        , emit: log_final
    tuple val(meta), path('*Log.out')                              , emit: log_out
    tuple val(meta), path('*Log.progress.out')                     , emit: log_progress
    path  "versions.yml"                                           , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.SJ.out.tab')            , optional:true, emit: sjdb

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sjdbFileChrStartEnd = sjdb ? "--sjdbFileChrStartEnd $sjdb" : ''
    def seq_center = meta.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$meta.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def out_sam_type = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    // def read_pair = params.protocol.contains("chromium") ? "${reads[1]} ${reads[0]}" : "${reads[0]} ${reads[1]}" -- commented out to be removed is it is not being used

    // default values max percentile for UMI count 0.99 and max to min ratio for UMI count 10 taken from STARsolo usage
    def cell_filter = meta.expected_cells ? "--soloCellFilter CellRanger2.2 $meta.expected_cells 0.99 10" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    if [[ $whitelist == *.gz ]]; then
        gzip -cdf $whitelist > whitelist.uncompressed.txt
    else
        cp $whitelist whitelist.uncompressed.txt
    fi

    STAR \\
        --genomeDir $index \\
        --readFilesIn ${(reverse ?: []).join( "," )} ${forward.join( "," )} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --sjdbGTFfile $gtf \\
        --soloCBwhitelist whitelist.uncompressed.txt \\
        --soloType $protocol \\
        --soloFeatures $star_feature \\
        --limitBAMsortRAM ${task.memory.toBytes()} \\
        $sjdbFileChrStartEnd \\
        $out_sam_type \\
        $seq_center \\
        $cell_filter \\
        $args \\

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    if [ -d ${prefix}.Solo.out ]; then
        # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi

    if [ -d ${prefix}.Solo.out/Velocyto/raw ]; then
        mv ${prefix}.Solo.out/Velocyto/raw ${prefix}.Solo.out/Velocyto/velocyto_raw
    fi

    if [ -d ${prefix}.Solo.out/Velocyto/filtered ]; then
        mv ${prefix}.Solo.out/Velocyto/filtered ${prefix}.Solo.out/Velocyto/velocyto_filtered
    fi

    if [ -d ${prefix}.Solo.out/SJ/raw ]; then 
        mv ${prefix}.Solo.out/SJ/raw ${prefix}.Solo.out/SJ/sj_raw
        gzip -f ${prefix}.Solo.out/SJ/sj_raw/features.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
