process createFaidxProcess {
    /*
    ==========================================================================
    CREATE FAIDX PROCESS
    ==========================================================================
    Creates FASTA index (.fai) for the reference genome
    Required by samtools for efficient random access

    Input:
        - fasta: Reference genome FASTA file

    Output:
        - fai: FASTA index file
    ==========================================================================
    */

    label 'process_single'
    container 'sageloom/gatk_minimal:latest'

    input:
        path fasta

    output:
        path("${fasta.baseName}.fai"), emit: fai
        path("${fasta.baseName}.amb"), emit: bwa_amb
        path("${fasta.baseName}.ann"), emit: bwa_ann
        path("${fasta.baseName}.bwt"), emit: bwa_bwt
        path("${fasta.baseName}.pac"), emit: bwa_pac
        path("${fasta.baseName}.sa"), emit: bwa_sa

    script:
        """
        # Create FASTA index (.fai) if it doesn't exist
        if [ ! -f ${fasta}.fai ]; then
            samtools faidx ${fasta}
        else
            cp ${fasta}.fai ${fasta}.fai.tmp
            mv ${fasta}.fai.tmp ${fasta}.fai
        fi

        # Create sequence index files for BWA if they don't exist
        if [ ! -f ${fasta}.bwt ]; then
            bwa index ${fasta}
        fi
        """
}
