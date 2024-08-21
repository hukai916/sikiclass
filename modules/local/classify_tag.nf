process CLASSIFY_TAG {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    path tag_fq

    output:
    tuple val(meta), path("01e_tmp_fasta/*.fasta"), emit: fasta
    tuple val(meta), path("01f_tmp_bam/*.bam"), emit: bam
    tuple val(meta), path("01a_no_tag/*.fastq.gz"), emit: no_tag
    tuple val(meta), path("01b_single_tag/*.fastq.gz"), emit: single_tag
    tuple val(meta), path("01c_multiple_tag/*.fastq.gz"), emit: multiple_tag
    tuple val(meta), path("01d_any_tag/*.fastq.gz"), emit: any_tag

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    # Step1: convert fastq to fasta
    mkdir 01e_tmp_fasta
    seqtk seq -a $reads > 01e_tmp_fasta/${prefix}.fasta

    # Step2: map tag to each individual read using bwa
    mkdir 01f_tmp_bam
    split_bwa.py 01e_tmp_fasta/${prefix}.fasta $tag_fq 01f_tmp_bam/${prefix}.bam

    # Step3: classify according to tag occurrence
    mkdir 01a_no_tag 01b_single_tag 01c_multiple_tag 01d_any_tag
    parse_bam_tag.py 01f_tmp_bam/${prefix}.bam $reads \
        01a_no_tag/${prefix}.fastq.gz \
        01b_single_tag/${prefix}.fastq.gz \
        01c_multiple_tag/${prefix}.fastq.gz \
        01d_any_tag/${prefix}.fastq.gz
    """

}
