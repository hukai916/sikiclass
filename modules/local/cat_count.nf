process CAT_COUNT {
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    path tsv
    val outdir
    val fq_class

    output:
    tuple val(fq_class), path("*/*.tsv"), emit: cat_count

    when:
    task.ext.when == null || task.ext.when

    // 00_stat 
        // 00a_fq_total, 00b_fq_no_tag_indel, 00c_fq_no_tag_del, 00d_fq_no_tag_insertion, 00e_fq_single_tag, 00f_fq_multiple_tag, 00g_fq_single_tag_precise, 00h_fq_single_tag_any_del, 00i_fq_single_tag_5del, 00j_fq_single_tag_3del, 00k_fq_single_tag_any_indel 000_master_table
    """
    mkdir $outdir
    cat *.tsv > ${outdir}/all_sample_${fq_class}.tsv

    """
}
