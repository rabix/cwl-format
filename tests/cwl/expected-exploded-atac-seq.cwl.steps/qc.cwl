#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: 'ATAC-seq 01 QC - reads: SE'

requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  default_adapters_file:
    doc: Adapters file
    type: File
  input_fastq_files:
    doc: Input fastq files
    type: File[]
  nthreads:
    doc: Number of threads.
    type: int

outputs:
  output_count_raw_reads:
    type: File[]
    outputSource: count_raw_reads/output_read_count
  output_custom_adapters:
    type: File[]
    outputSource: overrepresented_sequence_extract/output_custom_adapters
  output_diff_counts:
    type: File[]
    outputSource: compare_read_counts/result
  output_fastqc_data_files:
    doc: FastQC data files
    type: File[]
    outputSource: extract_fastqc_data/output_fastqc_data_file
  output_fastqc_report_files:
    doc: FastQC reports in zip format
    type: File[]
    outputSource: fastqc/output_qc_report_file

steps:
  compare_read_counts:
    in:
      file1: count_raw_reads/output_read_count
      file2: count_fastqc_reads/output_fastqc_read_count
    scatter:
    - file1
    - file2
    scatterMethod: dotproduct
    run: qc.cwl.steps/compare_read_counts.cwl
    out:
    - result
  count_fastqc_reads:
    in:
      input_basename: extract_basename/output_basename
      input_fastqc_data: extract_fastqc_data/output_fastqc_data_file
    scatter:
    - input_fastqc_data
    - input_basename
    scatterMethod: dotproduct
    run: qc.cwl.steps/count_fastqc_reads.cwl
    out:
    - output_fastqc_read_count
  count_raw_reads:
    in:
      input_basename: extract_basename/output_basename
      input_fastq_file: input_fastq_files
    scatter:
    - input_fastq_file
    - input_basename
    scatterMethod: dotproduct
    run: qc.cwl.steps/count_raw_reads.cwl
    out:
    - output_read_count
  extract_basename:
    in:
      input_file: input_fastq_files
    scatter: input_file
    run: qc.cwl.steps/extract_basename.cwl
    out:
    - output_basename
  extract_fastqc_data:
    in:
      input_basename: extract_basename/output_basename
      input_qc_report_file: fastqc/output_qc_report_file
    scatter:
    - input_qc_report_file
    - input_basename
    scatterMethod: dotproduct
    run: qc.cwl.steps/extract_fastqc_data.cwl
    out:
    - output_fastqc_data_file
  fastqc:
    in:
      input_fastq_file: input_fastq_files
      threads: nthreads
    scatter: input_fastq_file
    run: qc.cwl.steps/fastqc.cwl
    out:
    - output_qc_report_file
  overrepresented_sequence_extract:
    in:
      default_adapters_file: default_adapters_file
      input_basename: extract_basename/output_basename
      input_fastqc_data: extract_fastqc_data/output_fastqc_data_file
    scatter:
    - input_fastqc_data
    - input_basename
    scatterMethod: dotproduct
    run: qc.cwl.steps/overrepresented_sequence_extract.cwl
    out:
    - output_custom_adapters
