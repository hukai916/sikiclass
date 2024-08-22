process STAT_NO_TAG {
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    path tsv
    val pam 

    output:
    path "no_tag_indel_size_location_to_pam/*.tsv", emit: table_indel_size_location_to_pam
    path "no_tag_indel_size_distribution.tsv", emit: table_indel_size_distribution

    when:
    task.ext.when == null || task.ext.when

    // 00_stat/
    //  no_tag_indel_size_location_to_pam/*.tsv
    //  no_tag_indel_size_distribution.tsv
    """

    # Step1: make indel size and location to pam table
    mkdir -p no_tag_indel_size_location_to_pam
    for x in *.tsv; do
        if [ -e "\$x" ]; then
            outtsv=no_tag_indel_size_location_to_pam/\$x
            parse_indel_size_location.py \$x $pam \$outtsv
        fi
    done

    # Step2: make indel size distribution table 
    for x in *.tsv; do
        if [ -e "\$x" ]; then
            filename=\$(basename \$x .tsv)
            outtsv=no_tag_indel_size_distribution/\$x
            parse_indel_size_distribution.py \$x \$outtsv \${filename}
        fi
    done
    echo "sample\tdeletion_max\tdeletion_min\tdeletion_mean\tinsertion_max\tinsertion_min\tinsertion_mean" > no_tag_indel_size_distribution.tsv
    cat no_tag_indel_size_distribution/*.tsv >> no_tag_indel_size_distribution.tsv

    """
}
