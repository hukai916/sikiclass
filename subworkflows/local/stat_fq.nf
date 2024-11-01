include { COUNT_FQ } from '../../modules/local/count_fq'
include { CAT_COUNT } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ1 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT1 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ2 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT2 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ3 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT3 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ3A } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT3A } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ4 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT4 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ5 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT5 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ6 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT6 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ7 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT7 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ8 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT8 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ9 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT9 } from '../../modules/local/cat_count'
include { COUNT_FQ as COUNT_FQ10 } from '../../modules/local/count_fq'
include { CAT_COUNT as CAT_COUNT10 } from '../../modules/local/cat_count'
include { GET_TABLE_RATIO       } from '../../modules/local/get_table_ratio'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow STAT_FQ {

    take:
    fq_total
    fq_no_tag_indel
    fq_no_tag_deletion_only
    fq_no_tag_insertion_only
    fq_complex
    fq_single_tag
    fq_multiple_tag
    fq_any_tag
    fq_single_tag_precise_tag
    fq_single_tag_5indel
    fq_single_tag_3indel
    fq_single_tag_anyindel
    outfile

    main:
    ch_versions = Channel.empty()

    // Count fastq read numbers
    COUNT_FQ ( fq_total, "00a_fq_total" )
    CAT_COUNT ( COUNT_FQ.out.fq_c.collect(), "00a_fq_total", "fq_total" )
    COUNT_FQ1 ( fq_no_tag_indel, "00b_fq_no_tag_indel" )
    CAT_COUNT1 ( COUNT_FQ1.out.fq_c.collect(), "00b_fq_no_tag_indel", "fq_no_tag_indel" )
    COUNT_FQ2 ( fq_no_tag_deletion_only, "00c_fq_no_tag_deletion_only" )
    CAT_COUNT2 ( COUNT_FQ2.out.fq_c.collect(), "00c_fq_no_tag_deletion_only", "fq_no_tag_deletion_only" )
    COUNT_FQ3 ( fq_no_tag_insertion_only, "00d_fq_no_tag_insertion_only" )
    CAT_COUNT3 ( COUNT_FQ3.out.fq_c.collect(), "00d_fq_no_tag_insertion_only", "fq_no_tag_insertion_only" )
    COUNT_FQ3A ( fq_complex, "00e_fq_no_tag_complex" )
    CAT_COUNT3A ( COUNT_FQ3A.out.fq_c.collect(), "00e_fq_no_tag_complex", "fq_no_tag_complex" )
    COUNT_FQ4 ( fq_single_tag, "00f_fq_single_tag" )
    CAT_COUNT4 ( COUNT_FQ4.out.fq_c.collect(), "00f_fq_single_tag", "fq_single_tag" )
    COUNT_FQ5 ( fq_multiple_tag, "00g_fq_multiple_tag" )
    CAT_COUNT5 ( COUNT_FQ5.out.fq_c.collect(), "00g_fq_multiple_tag", "fq_multiple_tag" )
    COUNT_FQ6 ( fq_any_tag, "00h_fq_any_tag" )
    CAT_COUNT6 ( COUNT_FQ6.out.fq_c.collect(), "00h_fq_any_tag", "fq_any_tag" )
    COUNT_FQ7 ( fq_single_tag_precise_tag, "00i_fq_single_tag_precise" )
    CAT_COUNT7 ( COUNT_FQ7.out.fq_c.collect(), "00i_fq_single_tag_precise", "fq_single_tag_precise" )
    COUNT_FQ8 ( fq_single_tag_5indel, "00j_fq_single_tag_5del" )
    CAT_COUNT8 ( COUNT_FQ8.out.fq_c.collect(), "00j_fq_single_tag_5del", "fq_single_tag_5del" )
    COUNT_FQ9 ( fq_single_tag_3indel, "00k_fq_single_tag_3del" )
    CAT_COUNT9 ( COUNT_FQ9.out.fq_c.collect(), "00k_fq_single_tag_3del", "fq_single_tag_3del" )
    COUNT_FQ10 ( fq_single_tag_anyindel, "00l_fq_single_tag_any_indel" )
    CAT_COUNT10 ( COUNT_FQ10.out.fq_c.collect(), "00l_fq_single_tag_any_indel", "fq_single_tag_any_indel" )

    //  Make master tables 
        // 00_stat/fq_class_ratios.tsv
    GET_TABLE_RATIO (
        CAT_COUNT.out.cat_count,
        CAT_COUNT1.out.cat_count,
        CAT_COUNT2.out.cat_count,
        CAT_COUNT3.out.cat_count,
        CAT_COUNT3A.out.cat_count,
        CAT_COUNT4.out.cat_count,
        CAT_COUNT5.out.cat_count,
        CAT_COUNT6.out.cat_count,
        CAT_COUNT7.out.cat_count,
        CAT_COUNT8.out.cat_count,
        CAT_COUNT9.out.cat_count,
        CAT_COUNT10.out.cat_count,
        outfile
    )

}
