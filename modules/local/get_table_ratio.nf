process GET_TABLE_RATIO {
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    tuple val(total), path(tsv_total)
    tuple val(no_tag_indel), path(tsv_no_tag_indel)
    tuple val(no_tag_deletion), path(tsv_no_tag_deletion)
    tuple val(no_tag_insertion), path(tsv_no_tag_insertion)
    tuple val(no_tag_complex), path(tsv_no_tag_complex)
    tuple val(single_tag), path(tsv_single_tag)
    tuple val(multiple_tag), path(tsv_multiple_tag)
    tuple val(any_tag), path(tsv_any_tag)
    tuple val(single_tag_precise_tag), path(tsv_single_tag_precise_tag)
    tuple val(single_tag_5indel), path(tsv_single_tag_5indel)
    tuple val(single_tag_3indel), path(tsv_single_tag_3indel)
    tuple val(single_tag_anyindel), path(tsv_single_tag_anyindel)
    val outfile

    output:
    path "*.tsv", emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    get_table_ratio.py $outfile \
    $tsv_total $total \
    $tsv_no_tag_indel $no_tag_indel \
    $tsv_no_tag_deletion $no_tag_deletion \
    $tsv_no_tag_insertion $no_tag_insertion \
    $tsv_no_tag_complex $no_tag_complex \
    $tsv_single_tag $single_tag \
    $tsv_multiple_tag $multiple_tag \
    $tsv_any_tag $any_tag \
    $tsv_single_tag_precise_tag $single_tag_precise_tag \
    $tsv_single_tag_5indel $single_tag_5indel \
    $tsv_single_tag_3indel $single_tag_3indel \
    $tsv_single_tag_anyindel $single_tag_anyindel
    """

}
