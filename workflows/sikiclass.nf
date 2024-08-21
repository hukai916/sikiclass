/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CLASSIFY_TAG           } from '../modules/local/classify_tag'
include { CLASSIFY_SINGLE_TAG    } from '../modules/local/classify_single_tag'
include { CLASSIFY_NO_TAG        } from '../modules/local/classify_no_tag'
include { STAT_FQ                } from '../subworkflows/local/stat_fq'

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
    // 01_tag/
    //  01a_no_tag, 01b_single_tag, 01c_multiple_tag, 01d_any_tag, 01e_tmp_fasta, 01f_tmp_bam
    CLASSIFY_TAG ( 
        ch_samplesheet,
        tag_fq
    )

    // 
    // MODULE: classify reads with single tag with minimap2
    // 02_single_tag/
    //  02a_precise_tag, 02b_5indel, 02c_3indel, 02d_any_indel, 02e_tmp_bam, 02f_tmp_indel_pos
    CLASSIFY_SINGLE_TAG (
        CLASSIFY_TAG.out.single_tag,
        ref_with_tag,
        tag_start_ref_with_tag,
        tag_end_ref_with_tag,
        tag_flanking
    )

    // 
    // MODULE: classify reads without tag with minimap2
    // 03_no_tag/
    //  03a_indel, 03b_deletion, 03c_insertion, 03d_tmp_bam, 03e_tmp_indel_pos
    CLASSIFY_NO_TAG (
        CLASSIFY_TAG.out.no_tag,
        ref_wt,
        pam_start_ref_wt
    )


    // 
    // SUBWORKFLOW: stat fastq counts 
    // 00_stat/
    //  00a_fq_total, 00b_fq_no_tag_indel, 00c_fq_no_tag_del, 00d_fq_no_tag_insertion, 00e_fq_single_tag, 00f_fq_multiple_tag, 00g_fq_single_tag_precise, 00h_fq_single_tag_any_del, 00i_fq_single_tag_5del, 00j_fq_single_tag_3del, 00k_fq_single_tag_any_indel 000_master_table
    STAT_FQ (
        ch_samplesheet,
        CLASSIFY_NO_TAG.out.indel,
        CLASSIFY_NO_TAG.out.deletion,
        CLASSIFY_NO_TAG.out.insertion,
        CLASSIFY_TAG.out.single_tag,
        CLASSIFY_TAG.out.multiple_tag,
        CLASSIFY_TAG.out.any_tag,
        CLASSIFY_SINGLE_TAG.out.precise_tag,
        CLASSIFY_SINGLE_TAG.out.five_indel,
        CLASSIFY_SINGLE_TAG.out.three_indel,
        CLASSIFY_SINGLE_TAG.out.any_indel
    )


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
