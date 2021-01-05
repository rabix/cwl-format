#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: ATAC-seq - Quantification

requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  input_bam_files:
    type: File[]
  input_genome_sizes:
    type: File
  nthreads:
    type: int
    default: 1

outputs:
  bigwig_norm_files:
    doc: signal files of pileup reads in RPKM
    type: File[]
    outputSource: bamcoverage/output_bam_coverage
  bigwig_raw_files:
    doc: Raw reads bigWig (signal) files
    type: File[]
    outputSource: bdg2bw-raw/output_bigwig

steps:
  bamcoverage:
    in:
      bam: input_bam_files
      binSize:
        valueFrom: ${return 1}
      extendReads:
        valueFrom: ${return 200}
      normalizeUsing:
        valueFrom: RPKM
      numberOfProcessors: nthreads
      output_suffix:
        valueFrom: .rpkm.bw
    scatter: bam
    run: quant.cwl.steps/bamcoverage.cwl
    out:
    - output_bam_coverage
  bdg2bw-raw:
    in:
      bed_graph: bedsort_genomecov/bed_file_sorted
      genome_sizes: input_genome_sizes
      output_suffix:
        valueFrom: .raw.bw
    scatter: bed_graph
    run: quant.cwl.steps/bdg2bw-raw.cwl
    out:
    - output_bigwig
  bedsort_genomecov:
    in:
      bed_file: bedtools_genomecov/output_bedfile
    scatter: bed_file
    run: quant.cwl.steps/bedsort_genomecov.cwl
    out:
    - bed_file_sorted
  bedtools_genomecov:
    in:
      bg:
        valueFrom: ${return true}
      g: input_genome_sizes
      ibam: input_bam_files
    scatter: ibam
    run: quant.cwl.steps/bedtools_genomecov.cwl
    out:
    - output_bedfile
