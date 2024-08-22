process CLASSIFY_SINGLE_TAG {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    path ref
    val tag_start 
    val tag_end 
    val tag_flanking

    output:
      
    tuple val(meta), path("02a_precise_tag/*.fastq.gz"), emit: precise_tag
    tuple val(meta), path("02b_5indel/*.fastq.gz"), emit: five_indel
    tuple val(meta), path("02c_3indel/*.fastq.gz"), emit: three_indel
    tuple val(meta), path("02d_any_indel/*.fastq.gz"), emit: any_indel
    tuple val(meta), path("02e_tmp_bam/*.bam"), emit: bam
    path "02e_tmp_bam/*.bam", emit: bam_pure
    tuple val(meta), path("02f_tmp_indel_pos/*.tsv"), emit: indel_pos 
    path "02a_precise_tag/*.fastq.gz", emit: precise_tag_pure
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 02a_precise_tag, 02b_5indel, 02c_3indel, 02d_any_indel, 02e_tmp_bam, 02f_tmp_indel_pos

    """
    
    # Step1: map reads to ref_with_tag with minimap2
    mkdir 02e_tmp_bam
    minimap2 -ax map-pb $ref ${reads} | samtools sort --threads ${task.cpus} -o 02e_tmp_bam/${prefix}.bam

    # Step2: parse bam to extract indel positions
    mkdir 02f_tmp_indel_pos
    parse_bam_indel.py 02e_tmp_bam/${prefix}.bam $ref 02f_tmp_indel_pos/${prefix}.tsv

    # Step3: split reads into precise_insert, any_indel, 5_indel, 3_indel using #2 results
    mkdir 02a_precise_tag 02b_5indel 02c_3indel 02d_any_indel
    classify_single_tag.py 02f_tmp_indel_pos/${prefix}.tsv $reads \
    $tag_start \
    $tag_end \
    $tag_flanking \
    02a_precise_tag/${prefix}.fastq.gz \
    02b_5indel/${prefix}.fastq.gz \
    02c_3indel/${prefix}.fastq.gz \
    02d_any_indel/${prefix}.fastq.gz
    """

}
