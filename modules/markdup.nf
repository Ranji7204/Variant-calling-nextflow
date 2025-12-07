process markDuplicatesProcess {
    /*
    ==========================================================================
    MARK DUPLICATES PROCESS
    ==========================================================================
    Identifies and marks duplicate reads using Picard MarkDuplicates
    Creates index for the deduplicated BAM

    Input:
        - bam: Sorted BAM file from alignment
        - outdir: Output directory

    Output:
        - bam_dedup: BAM file with duplicates marked
        - metrics: Duplicate metrics file
    ==========================================================================
    */

    label 'process_high'
    container 'sageloom/gatk_minimal:latest'

    input:
        tuple val(sample_name), path(bam)
        path outdir

    output:
        tuple val(sample_name), path("${sample_name}_dedup.bam"), path("${sample_name}_dedup.bai"), emit: bam_dedup
        path("${sample_name}_dedup.metrics.txt"), emit: metrics

    script:
        """
        # Mark duplicates and create index
        picard MarkDuplicates \\
            I=${bam} \\
            O=${sample_name}_dedup.bam \\
            M=${sample_name}_dedup.metrics.txt \\
            CREATE_INDEX=true

        # Ensure BAI index exists with expected name
        if [ -f "${sample_name}_dedup.bai" ]; then
            echo "Index already in correct format"
        else
            # If Picard created .bam.bai, rename it
            if [ -f "${sample_name}_dedup.bam.bai" ]; then
                mv ${sample_name}_dedup.bam.bai ${sample_name}_dedup.bai
            fi
        fi

        # Copy to output directory
        cp ${sample_name}_dedup.bam ${outdir}/
        cp ${sample_name}_dedup.bai ${outdir}/
        cp ${sample_name}_dedup.metrics.txt ${outdir}/
        """
}
