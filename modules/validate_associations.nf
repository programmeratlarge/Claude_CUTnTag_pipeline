/*
 * Validate the user-provided association_csv against the FASTQ pairs TSV.
 * Enforces:
 *   - sample_id consistency
 *   - merge_group_id internal consistency (species/genome/antibody/condition)
 *   - control_group_id existence
 *   - peak_calling_mode validity
 */

process VALIDATE_ASSOCIATIONS {
    label 'process_low'
    publishDir "${params.outdir}/00_fastq_pairs", mode: params.publish_mode

    input:
    path pairs_tsv
    path association_csv
    val  allow_no_control
    val  genome

    output:
    path 'associations_validated.csv', emit: validated_csv
    path 'associations_validation.log', emit: log

    script:
    """
    validate_associations.py \\
        --pairs_tsv ${pairs_tsv} \\
        --csv       ${association_csv} \\
        --out_csv   associations_validated.csv \\
        --log       associations_validation.log \\
        --allow_no_control ${allow_no_control} \\
        --genome    '${genome}'
    """
}
