# GATK Mutect2 Nextflow Workflow

A production-ready Nextflow workflow for somatic variant calling using GATK Mutect2. This pipeline supports both tumor-only and tumor-normal analysis modes.

## Overview

This workflow implements the complete somatic variant calling pipeline:

```
FASTQ Files → Alignment (BWA) → MarkDuplicates → Mutect2 Calling → FilterMutectCalls → VCF
```

### Pipeline Steps

1. **Alignment** (`alignment.nf`)
   - BWA-MEM mapping to reference genome
   - Read group header assignment (sample name, platform)
   - Coordinate sorting and indexing
   - Processes all FASTQ files in parallel

2. **Mark Duplicates** (`markdup.nf`)
   - Picard MarkDuplicates to identify duplicate reads
   - Index generation for efficient access
   - Metrics file generation

3. **Reference Preparation**
   - **FAIDX** (`faidx.nf`): FASTA index creation for random access
   - **Sequence Dictionary** (`dictionary.nf`): GATK-required sequence dictionary

4. **Mutect2 Variant Calling** (`mutect2.nf`)
   - Somatic variant detection
   - Supports tumor-only mode (default)
   - Supports tumor-normal mode (with `--normal_sample` parameter)
   - Targets only regions in Baits.bed

5. **Variant Filtering** (`filter_mutect.nf`)
   - FilterMutectCalls for quality-based filtering
   - Removes low-confidence variants

## Data Structure

```
├── data/
│   ├── input/
│   │   ├── fastqs/
│   │   │   ├── chr13_R1.fastq.gz
│   │   │   ├── chr13_R2.fastq.gz
│   │   │   ├── chr17_R1.fastq.gz
│   │   │   ├── chr17_R2.fastq.gz
│   │   │   ├── chr3_R1.fastq.gz
│   │   │   ├── chr3_R2.fastq.gz
│   │   │   ├── chr4_R1.fastq.gz
│   │   │   ├── chr4_R2.fastq.gz
│   │   │   ├── chr7_R1.fastq.gz
│   │   │   └── chr7_R2.fastq.gz
│   │   └── Resources/
│   │       ├── Baits.bed
│   │       └── Baits.nochr.bed
│   ├── outputs/
│   │   └── (generated results)
│   └── References/
│       ├── hg19.p13.plusMT.no_alt_analysis_set.fa
│       ├── hg19.p13.plusMT.no_alt_analysis_set.fa.fai
│       ├── hg19.p13.plusMT.no_alt_analysis_set.dict
│       └── hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz
```

## Usage

### Basic Setup

1. Extract BWA index (if not already present):
```bash
cd training/data/References
tar -xzf hg19.p13.plusMT.no_alt_analysis_set.bwa_index.tar.gz
```

2. Navigate to workflow directory:
```bash
cd training/bk
```

### Run Modes

#### Tumor-Only Mode (Default)
```bash
nextflow run mutect2_workflow.nf \
    -c nextflow_mutect2.config \
    -params-file params_mutect2.yaml \
    -profile docker
```

#### Tumor-Normal Mode
```bash
nextflow run mutect2_workflow.nf \
    -c nextflow_mutect2.config \
    -params-file params_mutect2.yaml \
    -profile docker \
    --tumor_sample tumor_sample_name \
    --normal_sample normal_sample_name
```

#### With Custom Parameters
```bash
nextflow run mutect2_workflow.nf \
    -c nextflow_mutect2.config \
    -params-file params_mutect2.yaml \
    -profile docker \
    --tumor_sample my_tumor \
    --normal_sample my_normal \
    --output_dir ./custom_outputs
```

### Resume Previous Run
```bash
nextflow run mutect2_workflow.nf \
    -c nextflow_mutect2.config \
    -params-file params_mutect2.yaml \
    -profile docker \
    -resume
```

### View DAG
```bash
nextflow run mutect2_workflow.nf \
    -c nextflow_mutect2.config \
    -params-file params_mutect2.yaml \
    -profile docker \
    -preview
```

## Parameter Configuration

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tumor_sample` | test | Tumor sample name for Mutect2 |
| `normal_sample` | empty | Normal sample name (tumor-normal mode) |
| `reference_fasta` | hg19.p13.plusMT.no_alt_analysis_set.fa | Reference genome |
| `baits_bed` | Baits.bed | Target regions file |
| `docker_image` | sageloom/gatk_minimal:latest | Docker container image |
| `output_dir` | data/outputs | Output directory |

### Path Parameters

The workflow uses relative paths from the `bk/` directory:
- `base_dir`: Points to `../../data`
- All other paths are relative to `base_dir`

To use absolute paths, edit `nextflow_mutect2.config` or pass as parameters.

## Resource Allocation

The workflow is configured for systems with:
- **4 CPU cores** maximum
- **16 GB RAM** total

### Process Resources

| Process | CPUs | Memory | Label |
|---------|------|--------|-------|
| Alignment | 4 | 12 GB | process_high |
| MarkDuplicates | 4 | 12 GB | process_high |
| Mutect2 | 4 | 12 GB | process_high |
| FilterMutectCalls | 4 | 12 GB | process_high |
| FAIDX | 1 | 2 GB | process_single |
| Dictionary | 1 | 2 GB | process_single |

To adjust for different resource constraints, edit the `process` block in `nextflow_mutect2.config`.

## Output Files

### Per-Sample Outputs

For each input sample (e.g., `chr13`):

```
data/outputs/
├── chr13_dedup.bam              # Deduplicated BAM file
├── chr13_dedup.bai              # BAM index
├── chr13_dedup.metrics.txt       # Duplicate metrics
├── chr13.unfiltered.vcf.gz       # Raw Mutect2 calls
├── chr13.unfiltered.vcf.gz.tbi   # VCF index
├── chr13.filtered.vcf.gz         # Filtered final VCF
└── chr13.filtered.vcf.gz.tbi     # Filtered VCF index
```

### Workflow Reports

```
data/outputs/
├── execution_trace.txt           # Task execution trace
├── execution_report.html          # HTML execution report
├── execution_timeline.html        # Timeline visualization
└── execution_dag.svg              # Pipeline DAG
```

## Key Features

✅ **Parallel Processing**: All samples processed concurrently within resource constraints

✅ **Flexible Mode**: Switch between tumor-only and tumor-normal analysis with a parameter

✅ **Modular Design**: Each step is a separate process module for reusability

✅ **Docker Integration**: Pre-configured for containerized execution

✅ **Resource Management**: Automatic resource allocation based on available system resources

✅ **Comprehensive Reporting**: Execution traces, timelines, and DAGs generated automatically

✅ **Error Handling**: Workflow-level error handling and reporting

✅ **Caching**: Nextflow caching enabled for efficient re-runs

## Error Handling

If a process fails:

1. Check the error message in the terminal output
2. Examine task logs in `.nextflow/logs/`
3. Use `-resume` to restart from the failed point

Example:
```bash
nextflow run mutect2_workflow.nf -c nextflow_mutect2.config -profile docker -resume
```

## Dependencies

The `sageloom/gatk_minimal` Docker image must include:
- BWA (for alignment)
- samtools (for SAM/BAM operations)
- Picard tools (for MarkDuplicates, CreateSequenceDictionary)
- GATK (for Mutect2 and filtering)

## Directory Structure Reference

```
bk/                              # Main workflow directory
├── mutect2_workflow.nf           # Main workflow file
├── nextflow_mutect2.config       # Nextflow configuration
├── params_mutect2.yaml           # Parameters file
├── README_MUTECT2.md             # This file
└── modules/
    ├── alignment.nf              # BWA alignment process
    ├── markdup.nf                # MarkDuplicates process
    ├── faidx.nf                  # FASTA index process
    ├── dictionary.nf             # Sequence dictionary process
    ├── mutect2.nf                # Mutect2 calling process
    └── filter_mutect.nf          # FilterMutectCalls process
```

## Notes for Students

1. **Sample Naming**: FASTQ files are identified by their `_R1.fastq.gz` suffix. The basename becomes the sample name.

2. **Tumor-Normal Mode**: To use tumor-normal mode, you need BAM files from both tumor and normal samples. Ensure they are processed through the entire pipeline first.

3. **Multiple Samples**: The pipeline automatically discovers all FASTQ files in `data/input/fastqs/` and processes them in parallel.

4. **Resource Monitoring**: For systems with different resources, adjust the `process` blocks in `nextflow_mutect2.config`:
   ```groovy
   withLabel: process_high {
       cpus = 4      // Change this
       memory = '12 GB'  // And this
   }
   ```

5. **Clean Intermediate Files**: After successful completion, clean up with:
   ```bash
   nextflow clean -f
   ```

## Troubleshooting

### Common Issues

**Issue**: "Reference file not found"
- **Solution**: Ensure `data/References/` contains the FASTA file and BWA index

**Issue**: "R2 file not found for sample X"
- **Solution**: Check that paired FASTQ files follow the naming convention: `{sample_name}_R1.fastq.gz` and `{sample_name}_R2.fastq.gz`

**Issue**: Out of memory errors
- **Solution**: Reduce `task.cpus` or `memory` in `nextflow_mutect2.config`, or increase system memory

**Issue**: Docker image not found
- **Solution**: Pull the image: `docker pull sageloom/gatk_minimal:latest`

## References

- GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/articles/360035894731
- Mutect2 Documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360037593851
- Nextflow Documentation: https://www.nextflow.io/docs/latest/
- BWA: http://bio-bwa.sourceforge.net/

## License

This workflow is part of the Sageloom training repository.
