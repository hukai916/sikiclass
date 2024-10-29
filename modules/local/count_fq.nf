process COUNT_FQ {
    tag "$meta.id"
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(meta), path(reads)
    val outdir


    output:
    path "*/*.tsv", emit: fq_c

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 00_stat 
        // 00a_fq_total, 00b_fq_no_tag_indel, 00c_fq_no_tag_del, 00d_fq_no_tag_insertion, 00e_fq_single_tag, 00f_fq_multiple_tag, 00g_fq_single_tag_precise, 00h_fq_single_tag_any_del, 00i_fq_single_tag_5del, 00j_fq_single_tag_3del, 00k_fq_single_tag_any_indel 000_master_table
    """
    mkdir $outdir && count_fastq.py $reads $outdir/${prefix}.tsv $prefix

    """
}
