process GET_TABLE_RATIO {
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(total), path(csv_total)
    tuple val(no_tag_indel), path(csv_no_tag_indel)
    tuple val(no_tag_deletion), path(csv_no_tag_deletion)
    tuple val(no_tag_insertion), path(csv_no_tag_insertion)
    tuple val(single_tag), path(csv_single_tag)
    tuple val(multiple_tag), path(csv_multiple_tag)
    tuple val(any_tag), path(csv_any_tag)
    tuple val(single_tag_precise_tag), path(csv_single_tag_precise_tag)
    tuple val(single_tag_5indel), path(csv_single_tag_5indel)
    tuple val(single_tag_3indel), path(csv_single_tag_3indel)
    tuple val(single_tag_anyindel), path(csv_single_tag_anyindel)

    output:
    path "000_master_table/fq_class_ratios.tsv", emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir 000_master_table
    get_table_ratio.py 000_master_table/fq_class_ratios.tsv \
    $csv_total $total \
    $csv_no_tag_indel $no_tag_indel \
    $csv_no_tag_deletion $no_tag_deletion \
    $csv_no_tag_insertion $no_tag_insertion \
    $csv_single_tag $single_tag \
    $csv_multiple_tag $multiple_tag \
    $csv_any_tag $any_tag \
    $csv_single_tag_precise_tag $single_tag_precise_tag \
    $csv_single_tag_5indel $single_tag_5indel \
    $csv_single_tag_3indel $single_tag_3indel \
    $csv_single_tag_anyindel $single_tag_anyindel
    """

}
