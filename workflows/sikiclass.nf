/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CLASSIFY_TAG           } from '../modules/local/classify_tag'
include { CLASSIFY_SINGLE_TAG    } from '../modules/local/classify_single_tag'
include { CLASSIFY_NO_TAG        } from '../modules/local/classify_no_tag'
include { CLASSIFY_NO_TAG_FILTER_CONTROL } from '../modules/local/classify_no_tag_filter_control'
include { GET_CONTROL_INDEL      } from '../modules/local/get_control_indel'
include { STAT_FQ                } from '../subworkflows/local/stat_fq'
include { STAT_FQ as STAT_FQ_INDEL_FILTER_CONTROL } from '../subworkflows/local/stat_fq'
include { STAT_NO_TAG            } from '../modules/local/stat_no_tag'
include { STAT_NO_TAG as STAT_NO_TAG_FILTER_CONTROL } from '../modules/local/stat_no_tag'
include { STAT_SNP               } from '../modules/local/stat_snp'

// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sikiclass_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SIKICLASS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    tag_fq
    ref_with_tag
    ref_wt
    tag_start_ref_with_tag
    tag_end_ref_with_tag
    tag_flanking
    pam_start_ref_wt
    snp_pos 
    snp_wt
    snp_mut
    indel_range_to_scan_no_tag
    indel_range_to_scan_single_tag
    classify_no_tag_filter_control_indel

    main:

    ch_versions = Channel.empty()
    // Channel.fromPath( params.tag_fq ).set { tag_fq }
    // ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    // 00_fastqc/
    // FASTQC ( ch_samplesheet )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // 
    // MODULE: split reads according to tag occurence with BWA
    // 01_classify_tag/
    //  01a_no_tag, 01b_single_tag, 01c_multiple_tag, 01d_any_tag, 01e_tmp_fasta, 01f_tmp_bam
    CLASSIFY_TAG ( 
        ch_samplesheet,
        tag_fq
    )

    // 
    // MODULE: classify reads with single tag with minimap2
    // 02_classify_single_tag/
    //  02a_precise_tag, 02b_5indel, 02c_3indel, 02d_any_indel, 02e_tmp_bam, 02f_tmp_indel_pos
    CLASSIFY_SINGLE_TAG (
        CLASSIFY_TAG.out.single_tag,
        ref_with_tag,
        tag_start_ref_with_tag,
        tag_end_ref_with_tag,
        tag_flanking,
        indel_range_to_scan_single_tag
    )

    // 
    // MODULE: classify reads without tag with minimap2
    // 03_classify_no_tag/
    //  03a_indel, 03b_deletion_only, 03c_insertion_only, 03d_complex, 03e_tmp_bam, 03f_tmp_indel_pos
    CLASSIFY_NO_TAG (
        CLASSIFY_TAG.out.no_tag,
        ref_wt,
        pam_start_ref_wt,
        indel_range_to_scan_no_tag
    )

    // 
    // SUBWORKFLOW: stat fastq counts 
    // 00_stat/
    //  fq_class_ratios.tsv
    STAT_FQ (
        ch_samplesheet,
        CLASSIFY_NO_TAG.out.indel,
        CLASSIFY_NO_TAG.out.deletion_only,
        CLASSIFY_NO_TAG.out.insertion_only,
        CLASSIFY_NO_TAG.out.complex,
        CLASSIFY_TAG.out.single_tag,
        CLASSIFY_TAG.out.multiple_tag,
        CLASSIFY_TAG.out.any_tag,
        CLASSIFY_SINGLE_TAG.out.precise_tag,
        CLASSIFY_SINGLE_TAG.out.five_indel,
        CLASSIFY_SINGLE_TAG.out.three_indel,
        CLASSIFY_SINGLE_TAG.out.any_indel,
        "fq_class_ratios.tsv"
    )

    // 
    // MODULE: stat no tag reads regarding indel distribution 
    // 00_stat/
    //  no_tag_reads_indel_info/no_tag_indel_size_location_to_pam
    //  no_tag_reads_indel_info/no_tag_indel_size_distribution.tsv
    STAT_NO_TAG (
        CLASSIFY_NO_TAG.out.indel_pos_pure.collect(),
        pam_start_ref_wt,
        "no_tag_reads_indel_info"
    )

    // for "no tag" reads, if configured to filter out natural INDELS inferred from uninjected samples 
    if (classify_no_tag_filter_control_indel != null) {
        // 
        // MODULE: obtain natural INDELs from control samples
        // 03_classify_no_tag/
            // control_indels.tsv
        GET_CONTROL_INDEL ( CLASSIFY_NO_TAG.out.indel_pos_pure.collect(), classify_no_tag_filter_control_indel )

        // 
        // MODULE: classify reads without tag with minimap2 after excluding control indels, allow tuning minimum indel freq in control samples to be considered as real natural indels.
        // 03_classify_no_tag/
        //  03a_indel_filter_control, 03b_deletion_only_filter_control, 03c_insertion_only_filter_control, 03d_complex_filter_control, 03f_tmp_indel_pos_filter_control
        CLASSIFY_NO_TAG_FILTER_CONTROL (
            CLASSIFY_TAG.out.no_tag,
            ref_wt,
            pam_start_ref_wt,
            GET_CONTROL_INDEL.out.tsv
        )

        // 
        // MODULE: stat no tag reads regarding indel distribution 
        // 00_stat/
        //  no_tag_reads_indel_info_filter_control/no_tag_indel_size_location_to_pam_filter_control
        //  no_tag_reads_indel_info_filter_control/no_tag_indel_size_distribution_filter_control.tsv
        STAT_NO_TAG_FILTER_CONTROL (
            CLASSIFY_NO_TAG_FILTER_CONTROL.out.indel_pos_pure.collect(),
            pam_start_ref_wt,
            "no_tag_reads_indel_info_filter_control"
        )

        // 
        // SUBWORKFLOW: stat fastq counts 
        // 00_stat/
        //  fq_class_ratios_indel_filter_control.tsv
        STAT_FQ_INDEL_FILTER_CONTROL (
            ch_samplesheet,
            CLASSIFY_NO_TAG_FILTER_CONTROL.out.indel,
            CLASSIFY_NO_TAG_FILTER_CONTROL.out.deletion_only,
            CLASSIFY_NO_TAG_FILTER_CONTROL.out.insertion_only,
            CLASSIFY_NO_TAG_FILTER_CONTROL.out.complex,
            CLASSIFY_TAG.out.single_tag,
            CLASSIFY_TAG.out.multiple_tag,
            CLASSIFY_TAG.out.any_tag,
            CLASSIFY_SINGLE_TAG.out.precise_tag,
            CLASSIFY_SINGLE_TAG.out.five_indel,
            CLASSIFY_SINGLE_TAG.out.three_indel,
            CLASSIFY_SINGLE_TAG.out.any_indel,
            "fq_class_ratios_indel_filter_control.tsv"
        )
    }


    // 
    // MODULE: gather natural INDELs 
    // 03_classify_no_tag/
    //  03


    // 
    // MODULE: classify reads without tag with minimap2
    // 03_classify_no_tag/
    //  03a_indel, 03b_deletion, 03c_insertion, 03d_tmp_bam, 03e_tmp_indel_pos
    // CLASSIFY_NO_TAG (
    //     CLASSIFY_TAG.out.no_tag,
    //     ref_wt,
    //     pam_start_ref_wt
    // )

    if (snp_pos != null) {
        log.info "params.snp_pos supplied!"
        // MODULE: split and count SNP occurrences for precise_tag reads
        // 00_stat/000_master_table/precise_tag_snp_fraction.tsv
        // 02_classify_single_tag/02a_precise_tag_snp_wt/*.fastq.gz
        // 02_classify_single_tag/02a_precise_tag_snp_mut/*.fastq.gz
        // 02_classify_single_tag/02a_precise_tag_snp_other/*.fastq.gz
        STAT_SNP (
            CLASSIFY_SINGLE_TAG.out.precise_tag_pure.collect(),
            CLASSIFY_SINGLE_TAG.out.bam_pure.collect(),
            snp_pos,
            snp_wt,
            snp_mut
        )
        
    } else {
        log.info "params.snp_pos not supplied, end!"
    }

    // 
    // MODULE: stat fastq counts 
    // 00_stat/
    //  00a_fq_total, 00b_fq_no_tag_indel, 00c_fq_no_tag_del, 00d_fq_no_tag_insertion, 00e_fq_single_tag, 00f_fq_multiple_tag, 00g_fq_single_tag_precise, 00h_fq_single_tag_any_del, 00i_fq_single_tag_3del, 00j_fq_single_tag_5del, 000_master_table
    // STAT_FQ (

    // )


    //
    // Collate and save software versions
    //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
    //     .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    // ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    // ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    // summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )

    emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
