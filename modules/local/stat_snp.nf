process STAT_SNP {
    label 'process_medium'
    container 'hukai916/sikiclass:0.2'

    input:
    path fastq 
    path bam
    val snp_pos 
    val snp_wt
    val snp_mut

    output:
    path "00_stat/precise_tag_snp_fraction.tsv", emit: table_precise_tag_snp_fraction
    path "02_classify_single_tag/02a_precise_tag_snp_wt/*.fastq.gz", emit: snp_wt 
    path "02_classify_single_tag/02a_precise_tag_snp_mut/*.fastq.gz", emit: snp_mut 
    path "02_classify_single_tag/02a_precise_tag_snp_other/*.fastq.gz", emit: snp_other 

    when:
    task.ext.when == null || task.ext.when

    // 00_stat/precise_tag_snp_fraction.tsv
    // 02_classify_single_tag/02a_precise_tag_snp_wt/*.fastq.gz
    // 02_classify_single_tag/02a_precise_tag_snp_mut/*.fastq.gz
    // 02_classify_single_tag/02a_precise_tag_snp_other/*.fastq.gz
    """

    # Step1: make indel size and location to pam table
    mkdir -p 00_stat
    mkdir -p 02_classify_single_tag/02a_precise_tag_snp_wt
    mkdir -p 02_classify_single_tag/02a_precise_tag_snp_mut
    mkdir -p 02_classify_single_tag/02a_precise_tag_snp_other
    mkdir tmp

    for x in *.fastq.gz; do
        filename=\$(basename \$x .fastq.gz)
        bam=\${filename}.bam
        out_wt=02_classify_single_tag/02a_precise_tag_snp_wt/\${filename}.fastq.gz
        out_mut=02_classify_single_tag/02a_precise_tag_snp_mut/\${filename}.fastq.gz
        out_other=02_classify_single_tag/02a_precise_tag_snp_other/\${filename}.fastq.gz
        tsv=tmp/\${filename}.tsv
        parse_snp.py \$bam \$x $snp_pos $snp_wt $snp_mut \$out_wt \$out_mut \$out_other \$tsv \${filename}
    done

    echo "sample\tprecise_insert_fragments\twt_SNP_fragments\twt_SNP_ratio\tmut_SNP_fragments\tmut_SNP_ratio\tother_SNP_fragments\tother_SNP_ratio" > 00_stat/precise_tag_snp_fraction.tsv
    cat tmp/*.tsv >> 00_stat/precise_tag_snp_fraction.tsv

    """
}
