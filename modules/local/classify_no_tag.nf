process CLASSIFY_NO_TAG {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    path ref
    val pam_start

    output:
      
    tuple val(meta), path("03a_indel/*.fastq.gz"), emit: indel
    tuple val(meta), path("03b_deletion/*.fastq.gz"), emit: deletion
    tuple val(meta), path("03c_insertion/*.fastq.gz"), emit: insertion
    tuple val(meta), path("03d_tmp_bam/*.bam"), emit: bam
    tuple val(meta), path("03e_tmp_indel_pos/*.tsv"), emit: indel_pos
    path "03e_tmp_indel_pos/*.tsv", emit: indel_pos_pure
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 03a_indel, 03b_deletion, 03c_insertion, 03d_tmp_bam, 03e_tmp_indel_pos
    """
    
    # Step1: map reads to ref_wt with minimap2
    mkdir 03d_tmp_bam
    minimap2 -ax map-pb $ref ${reads} | samtools sort --threads ${task.cpus} -o 03d_tmp_bam/${prefix}.bam

    # Step2: parse bam to extract indel positions
    mkdir 03e_tmp_indel_pos
    parse_bam_indel.py 03d_tmp_bam/${prefix}.bam $ref 03e_tmp_indel_pos/${prefix}.tsv

    # Step3: split reads into indel, deletion, and insertion using #2 results
    mkdir 03a_indel 03b_deletion 03c_insertion
    classify_no_tag.py 03e_tmp_indel_pos/${prefix}.tsv $reads \
    $pam_start \
    03a_indel/${prefix}.fastq.gz \
    03b_deletion/${prefix}.fastq.gz \
    03c_insertion/${prefix}.fastq.gz
    """
}
