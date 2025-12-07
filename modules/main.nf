#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ==============================================================================
// MUTECT2 GENOMICS WORKFLOW
// ==============================================================================
// Performs DNA variant calling using GATK Mutect2
// Supports tumor-only and tumor-normal modes
// Process: Alignment -> MarkDuplicates -> Dictionary Creation -> Variant Calling -> Filtering
// ==============================================================================

include { alignmentProcess } from './modules/alignment.nf'
include { markDuplicatesProcess } from './modules/markdup.nf'
include { createFaidxProcess } from './modules/faidx.nf'
include { createDictionaryProcess } from './modules/dictionary.nf'
include { mutect2CallingProcess } from './modules/mutect2.nf'
include { filterMutectCallsProcess } from './modules/filter_mutect.nf'

workflow {
    // ===========================================================================
    // INPUT VALIDATION
    // ===========================================================================

    // Validate required parameters
    if (!params.fastq_dir) {
        error("Missing required parameter: --fastq_dir (directory containing FASTQ files)")
    }
    if (!params.ref_fa) {
        error("Missing required parameter: --ref_fa (path to reference FASTA file)")
    }
    if (!params.baits) {
        error("Missing required parameter: --baits (path to target regions BED file)")
    }
    if (!params.outdir) {
        error("Missing required parameter: --outdir (output directory)")
    }

    // Convert to file objects and validate
    fastq_dir = file(params.fastq_dir)
    ref_fasta_file = file(params.ref_fa)
    baits_bed_file = file(params.baits)
    outdir = file(params.outdir)

    // Keep absolute path strings for passing to processes in Docker
    ref_fasta_path = params.ref_fa
    baits_bed_path = params.baits

    if (!fastq_dir.exists()) {
        error("FASTQ directory does not exist: ${params.fastq_dir}")
    }
    if (!ref_fasta_file.exists()) {
        error("Reference FASTA file does not exist: ${params.ref_fa}")
    }
    if (!baits_bed_file.exists()) {
        error("Baits BED file does not exist: ${params.baits}")
    }

    // Create output directory if it doesn't exist
    outdir.mkdirs()

    // ===========================================================================
    // INPUT SETUP
    // ===========================================================================

    // Get all R1 fastq files from the provided directory
    fastq_files = Channel.fromPath("${fastq_dir}/*_R1.fastq.gz")
        .map { r1_file ->
            def sample_name = r1_file.name.replaceAll('_R1.fastq.gz', '')
            def r2_file = file("${fastq_dir}/${sample_name}_R2.fastq.gz")

            if (!r2_file.exists()) {
                error("R2 file not found for ${sample_name}: ${fastq_dir}/${sample_name}_R2.fastq.gz")
            }

            [sample_name, r1_file, r2_file]
        }

    main:
        // =======================
        // Reference detection and alignment
        // We expect prebuilt index files (.fai and .dict) to be present next to
        // the provided `--ref_fa`. If they are missing, fail early with a clear
        // message so the user can prepare them once (no repeated param).
        // =======================

        // Look for index files in the same directory as the reference FASTA
        // Use params.ref_fa string here to avoid calling methods on Nextflow `file()` wrapper
        ref_dir = new File(params.ref_fa).parentFile

        // Prefer exact .fai alongside the fasta (e.g. /path/ref.fa.fai)
        def fai_candidate = file("${params.ref_fa}.fai")
        def fai_file = fai_candidate.exists() ? fai_candidate : ref_dir.listFiles().find { it.name.endsWith('.fai') }
        if (!fai_file) {
            error("Missing FASTA index (.fai) next to reference FASTA: ${ref_fasta_file}. Please create the index and retry.")
        }

        // Find a sequence dictionary in the same directory
        def dict_file_local = ref_dir.listFiles().find { it.name.endsWith('.dict') }
        if (!dict_file_local) {
            error("Missing sequence dictionary (.dict) next to reference FASTA: ${ref_fasta_file}. Please create the dictionary and retry.")
        }

        // Broadcast the index files as a value channel so each alignment task gets them staged
        // Emit string paths to avoid file/path type issues
        dict_file = Channel.value(dict_file_local.toString())

        // Find BWA index files. Some environments create files like `ref.fa.pac` (basename includes .fa)
        // while others use `ref.pac`. Try both possibilities and accept whichever exists.
        def refName = new File(params.ref_fa).name
        def refBaseNoExt = refName.replaceAll(/(\.fa|\.fasta)(\.gz)?$/, '')
        def bwa_exts = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        def bwa_index_files = []
        def missing = []
        bwa_exts.each { ext ->
            def f1 = new File(ref_dir, refName + ext)
            def f2 = new File(ref_dir, refBaseNoExt + ext)
            if (f1.exists()) bwa_index_files << f1
            else if (f2.exists()) bwa_index_files << f2
            else missing << (refName + ext)
        }
        if (missing) {
            error("Missing BWA index files next to reference FASTA: ${missing}. Please run 'bwa index' on the reference or provide the index files alongside the FASTA.")
        }

        // Include the .fai and all BWA index files in a single channel so Nextflow stages them
        def all_index_files = ([fai_file] + bwa_index_files).collect { it.toString() }
        ref_index_ch = Channel.value(all_index_files)

        // Run alignment using the found index files
        alignmentProcess(fastq_files, ref_fasta_file, ref_index_ch)

        // Deduplication
        markDuplicatesProcess(alignmentProcess.out.bam, outdir)

        // Convert normal_sample to a channel value
        normal_sample = Channel.value(params.normal_sample ?: "")

        mutect2CallingProcess(
            markDuplicatesProcess.out.bam_dedup,
            ref_fasta_file,
            ref_index_ch,
            baits_bed_file,
            dict_file,
            normal_sample,
            outdir
        )

        // =======================
        // Filtering
        // =======================
        filterMutectCallsProcess(mutect2CallingProcess.out.vcf.map { sample_name, vcf, tbi, stats -> [sample_name, vcf, tbi, stats] }, ref_fasta_file, ref_index_ch, dict_file, outdir)

    emit:
        final_vcf = filterMutectCallsProcess.out.filtered_vcf
}

workflow.onComplete {
    println """
    =================================================================
    Pipeline execution summary
    =================================================================
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Results in  : ${params.outdir}
    Work dir    : ${workflow.workDir}
    Exit status : ${workflow.exitStatus}
    =================================================================
    """
}

workflow.onError {
    println """
    =================================================================
    Pipeline execution error
    =================================================================
    Error message: ${workflow.errorMessage}
    Error report : ${workflow.errorReport}
    =================================================================
    """
}
