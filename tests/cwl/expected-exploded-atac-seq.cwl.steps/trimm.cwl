#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: 'ATAC-seq 02 trimming - reads: SE'

requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  input_adapters_files:
    doc: Input adapters files
    type: File[]
  input_fastq_files:
    doc: Input fastq files
    type: File[]
  nthreads:
    doc: Number of threads
    type: int
    default: 1
  quality_score:
    type: string
    default: -phred33
  trimmomatic_jar_path:
    doc: Trimmomatic Java jar file
    type: string
    default: /usr/share/java/trimmomatic.jar
  trimmomatic_java_opts:
    doc: JVM arguments should be a quoted, space separated list
    type: string?

outputs:
  output_data_fastq_trimmed_files:
    doc: Trimmed fastq files
    type: File[]
    outputSource: trimmomatic/output_read1_trimmed_file
  output_trimmed_fastq_read_count:
    doc: Trimmed read counts of fastq files
    type: File[]
    outputSource: count_fastq_reads/output_read_count

steps:
  count_fastq_reads:
    in:
      input_basename: extract_basename/output_basename
      input_fastq_file: trimmomatic/output_read1_trimmed_file
    scatter:
    - input_fastq_file
    - input_basename
    scatterMethod: dotproduct
    run: trimm.cwl.steps/count_fastq_reads.cwl
    out:
    - output_read_count
  extract_basename:
    in:
      input_file: trimmomatic/output_read1_trimmed_file
    scatter: input_file
    run: trimm.cwl.steps/extract_basename.cwl
    out:
    - output_basename
  trimmomatic:
    in:
      end_mode:
        valueFrom: SE
      illuminaclip:
        valueFrom: 2:30:15
      input_adapters_file: input_adapters_files
      input_read1_fastq_file: input_fastq_files
      java_opts: trimmomatic_java_opts
      leading:
        valueFrom: ${return 3}
      minlen:
        valueFrom: ${return 15}
      nthreads: nthreads
      phred:
        valueFrom: '33'
      slidingwindow:
        valueFrom: 4:20
      trailing:
        valueFrom: ${return 3}
      trimmomatic_jar_path: trimmomatic_jar_path
    scatter:
    - input_read1_fastq_file
    - input_adapters_file
    scatterMethod: dotproduct
    run: trimm.cwl.steps/trimmomatic.cwl
    out:
    - output_read1_trimmed_file
