process filterMutectCallsProcess {
    /*
    ==========================================================================
    FILTER MUTECT CALLS PROCESS
    ==========================================================================
    Filters Mutect2 variant calls based on GATK recommended filters
    Removes low-confidence variants

    Input:
        - vcf: Unfiltered VCF from Mutect2
        - fasta: Reference genome FASTA
        - dict: Sequence dictionary

    Output:
        - filtered_vcf: Filtered VCF file
    ==========================================================================
    */

    label 'process_high'
    container 'sageloom/gatk_minimal:latest'

    input:
        tuple val(sample_name), path(vcf), path(tbi), path(stats)
        path fasta
        // Ensure fasta index is staged for GATK
        path index_files
        path dict
        path outdir

    output:
        tuple val(sample_name), path("${sample_name}.filtered.vcf.gz"), path("${sample_name}.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
        """
        # Filter Mutect2 calls
        gatk --java-options '-Xmx${task.memory.toGiga()}G' FilterMutectCalls \\
            -R ${fasta} \\
            -V ${vcf} \\
            -O ${sample_name}.filtered.vcf.gz

        # Verify index was created
        if [ ! -f "${sample_name}.filtered.vcf.gz.tbi" ]; then
            gatk IndexFeatureFile \\
                -I ${sample_name}.filtered.vcf.gz
        fi

        # Copy final outputs to output directory
        mkdir -p ${outdir}
        cp ${sample_name}.filtered.vcf.gz ${outdir}/
        cp ${sample_name}.filtered.vcf.gz.tbi ${outdir}/

        # Also copy the unfiltered VCF for reference
        cp ${vcf} ${outdir}/ || true
        cp ${tbi} ${outdir}/ || true
        """
}
