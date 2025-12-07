process alignmentProcess {
    /*
    ==========================================================================
    ALIGNMENT PROCESS
    ==========================================================================
    Performs read alignment using BWA MEM
    Outputs BAM file with index

    Input:
        - sample_name: Name of the sample (from fastq prefix)
        - fastq_r1: Forward reads (R1)
        - fastq_r2: Reverse reads (R2)
        - fasta: Reference genome FASTA
        - index_files: List of index files (staged into workdir)

    Output:
        - bam: Sorted and indexed BAM file
    ==========================================================================
    */

    label 'process_high'
    container 'sageloom/gatk_minimal:latest'

    input:
        tuple val(sample_name), path(fastq_r1), path(fastq_r2)
        path fasta
        // `index_files` is a list of index files (.fai, .amb, .ann, .bwt, .pac, .sa)
        path index_files

    output:
        tuple val(sample_name), path("${sample_name}.bam"), emit: bam

    script:
        """
        # Run BWA MEM alignment
        # -M: Mark shorter split hits as secondary (for compatibility with Picard)
        # -R: Read group header (includes sample name and platform)
        # -t: Number of threads

        bwa mem \\
            -M \\
            -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA" \\
            -t ${task.cpus} \\
            ${fasta} \\
            ${fastq_r1} ${fastq_r2} \\
            | samtools view -b -@ ${task.cpus} - \\
            | samtools sort -@ ${task.cpus} -o ${sample_name}.bam

        # Index the BAM file
        samtools index ${sample_name}.bam
        """
}
