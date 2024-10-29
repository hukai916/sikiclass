process CLASSIFY_NO_TAG {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    path ref
    val pam_start
    val indel_range_to_scan

    output:
      
    tuple val(meta), path("03a_indel/*.fastq.gz"), emit: indel
    tuple val(meta), path("03b_deletion_only/*.fastq.gz"), emit: deletion_only
    tuple val(meta), path("03c_insertion_only/*.fastq.gz"), emit: insertion_only 
    tuple val(meta), path("03d_complex/*.fastq.gz"), emit: complex 
    tuple val(meta), path("03e_tmp_bam/*.bam"), emit: bam
    tuple val(meta), path("03f_tmp_indel_pos/*.tsv"), emit: indel_pos
    path "03f_tmp_indel_pos/*.tsv", emit: indel_pos_pure
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    # Step1: map reads to ref_wt with minimap2
    mkdir 03e_tmp_bam && minimap2 -ax map-pb $ref ${reads} | samtools sort --threads ${task.cpus} -o 03e_tmp_bam/${prefix}.bam

    # Step2: parse bam to extract indel positions
    mkdir 03f_tmp_indel_pos && parse_bam_indel.py 03e_tmp_bam/${prefix}.bam $ref 03f_tmp_indel_pos/${prefix}.tem.tsv

    # Step2-1: filter indel positions by range
    filter_indel_by_range.py 03f_tmp_indel_pos/${prefix}.tem.tsv "$indel_range_to_scan" 03f_tmp_indel_pos/${prefix}.tsv

    # Step3: split reads into indel, deletion, and insertion using #2 results
    mkdir 03a_indel 03b_deletion_only 03c_insertion_only 03d_complex && \
    classify_no_tag.py 03f_tmp_indel_pos/${prefix}.tsv $reads \
    $pam_start \
    03a_indel/${prefix}.fastq.gz \
    03b_deletion_only/${prefix}.fastq.gz \
    03c_insertion_only/${prefix}.fastq.gz \
    03d_complex/${prefix}.fastq.gz
    """
}
