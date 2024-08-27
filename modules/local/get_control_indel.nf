process GET_CONTROL_INDEL {
    label 'process_medium'

    container 'hukai916/sikiclass:0.2'

    input:
    path reads
    val control_samples

    output:
    path "control_indels.tsv", emit: tsv
    
    when:
    task.ext.when == null || task.ext.when

    script:

    """
    get_control_indel.py "$control_samples" control_indels.tsv

    """
}
