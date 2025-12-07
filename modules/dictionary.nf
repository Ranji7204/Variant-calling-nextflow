process createDictionaryProcess {
    /*
    ==========================================================================
    CREATE SEQUENCE DICTIONARY PROCESS
    ==========================================================================
    Creates a sequence dictionary (.dict) from the reference FASTA
    Required by GATK for variant calling

    Input:
        - fasta: Reference genome FASTA file

    Output:
        - dict: Sequence dictionary file
    ==========================================================================
    */

    label 'process_single'
    container 'sageloom/gatk_minimal:latest'

    input:
        path fasta

    output:
        path("*.dict"), emit: dict

    script:
        // compute basename from the staged fasta path
        def dict_name = file(fasta).baseName + '.dict'

        """
        # Create sequence dictionary if it doesn't exist
        if [ ! -f ${dict_name} ]; then
            picard CreateSequenceDictionary \\
                REFERENCE=${fasta} \\
                OUTPUT=${dict_name}
        else
            # Copy it to ensure it's in the work directory
            cp ${dict_name} ${dict_name}.tmp
            mv ${dict_name}.tmp ${dict_name}
        fi
        """
}
