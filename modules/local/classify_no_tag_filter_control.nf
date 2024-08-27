process CLASSIFY_NO_TAG_FILTER_CONTROL {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    path ref
    val pam_start
    path control_indel

    output:
    tuple val(meta), path("03a_indel_filter_control/*.fastq.gz"), emit: indel
    tuple val(meta), path("03b_deletion_only_filter_control/*.fastq.gz"), emit: deletion_only
    tuple val(meta), path("03c_insertion_only_filter_control/*.fastq.gz"), emit: insertion_only 
    tuple val(meta), path("03d_complex_filter_control/*.fastq.gz"), emit: complex 
    // tuple val(meta), path("03e_tmp_bam/*.bam"), emit: bam
    tuple val(meta), path("03f_tmp_indel_pos_filter_control/*.tsv"), emit: indel_pos
    path "03f_tmp_indel_pos_filter_control/*.tsv", emit: indel_pos_pure
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control_indel_freq_cutoff = task.ext.control_indel_freq_cutoff ?: 1

    """
    
    # Step1: filter out control indels from raw indel tsv
    mkdir 03e_tmp_bam
    minimap2 -ax map-pb $ref ${reads} | samtools sort --threads ${task.cpus} -o 03e_tmp_bam/${prefix}.bam

    # Step2: parse bam to extract indel positions
    mkdir 03f_tmp_indel_pos
    parse_bam_indel.py 03e_tmp_bam/${prefix}.bam $ref 03f_tmp_indel_pos/${prefix}.tsv

    # Step3: filter out control INDELs
    mkdir 03f_tmp_indel_pos_filter_control
    filter_control_indel.py $control_indel 03f_tmp_indel_pos/${prefix}.tsv 03f_tmp_indel_pos_filter_control/${prefix}.tsv $prefix $control_indel_freq_cutoff

    # Step4: split reads into indel, deletion, and insertion using #3 results
    mkdir 03a_indel_filter_control 03b_deletion_only_filter_control 03c_insertion_only_filter_control 03d_complex_filter_control
    classify_no_tag.py 03f_tmp_indel_pos_filter_control/${prefix}.tsv $reads \
    $pam_start \
    03a_indel_filter_control/${prefix}.fastq.gz \
    03b_deletion_only_filter_control/${prefix}.fastq.gz \
    03c_insertion_only_filter_control/${prefix}.fastq.gz \
    03d_complex_filter_control/${prefix}.fastq.gz
    """
}
