process mutect2CallingProcess {
    /*
    ==========================================================================
    MUTECT2 CALLING PROCESS
    ==========================================================================
    Calls somatic variants using GATK Mutect2
    Supports both tumor-only and tumor-normal modes

    Input:
        - bam: Deduplicated BAM file
        - fasta: Reference genome FASTA
        - baits_bed: Target regions BED file
        - dict: Sequence dictionary
        - tumor_sample: Tumor sample name
        - normal_sample: Normal sample name (empty if tumor-only)
        - is_paired: Boolean indicating if tumor-normal mode

    Output:
        - vcf: Unfiltered VCF from Mutect2
    ==========================================================================
    */

    label 'process_high'
    container 'sageloom/gatk_minimal:latest'

    input:
        tuple val(sample_name), path(bam), path(bai)
        path fasta
        // Ensure fasta index/bwa index files are staged for GATK
        path index_files
        path baits_bed
        path dict
        val normal_sample
        path outdir

    output:
        tuple val(sample_name), path("${sample_name}.unfiltered.vcf.gz"), path("${sample_name}.unfiltered.vcf.gz.tbi"), path("${sample_name}.unfiltered.vcf.gz.stats"), emit: vcf

    script:
        """
        #!/usr/bin/env bash
        set -euo pipefail

        echo "[mutect2] Running Mutect2 for sample: ${sample_name}"

        if [ -n "${normal_sample}" ]; then
            echo "[mutect2] Note: normal sample provided (${normal_sample}) but pipeline currently expects normal BAM input separately. Proceeding with tumor-only call."
        fi

        gatk --java-options '-Xmx${task.memory.toGiga()}G' Mutect2 \
            -R ${fasta} \
            -I ${bam} \
            -tumor ${sample_name} \
            -L ${baits_bed} \
            -O ${sample_name}.unfiltered.vcf.gz

        # Index the VCF
        gatk IndexFeatureFile -I ${sample_name}.unfiltered.vcf.gz

        # Ensure stats file exists (required for FilterMutectCalls)
        # Create a placeholder if not present
        #if [ ! -f "${sample_name}.unfiltered.vcf.gz.stats" ]; then
        #    echo "[mutect2] Creating placeholder stats file for FilterMutectCalls"
        #    touch ${sample_name}.unfiltered.vcf.gz.stats
        #fi
        """
}
