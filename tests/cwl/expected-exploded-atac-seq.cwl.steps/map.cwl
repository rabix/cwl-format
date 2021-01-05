#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
doc: 'ATAC-seq 03 mapping - reads: SE'

requirements:
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  genome_ref_first_index_file:
    doc: |-
      Bowtie first index files for reference genome (e.g. *1.ebwt). The rest of the files should be in the same folder.
    type: File
    secondaryFiles:
    - ^^.2.ebwt
    - ^^.3.ebwt
    - ^^.4.ebwt
    - ^^.rev.1.ebwt
    - ^^.rev.2.ebwt
  genome_sizes_file:
    doc: Genome sizes tab-delimited file (used in samtools)
    type: File
  input_fastq_files:
    doc: Input fastq files
    type: File[]
  nthreads:
    type: int
    default: 1
  picard_jar_path:
    doc: Picard Java jar file
    type: string
    default: /usr/picard/picard.jar
  picard_java_opts:
    doc: |-
      JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
    type: string?

outputs:
  output_bowtie_log:
    doc: Bowtie log file.
    type: File[]
    outputSource: bowtie-se/output_bowtie_log
  output_data_sorted_dedup_bam_files:
    doc: BAM files without duplicate reads.
    type: File[]
    outputSource: index_dedup_bams/indexed_file
  output_data_sorted_dups_marked_bam_files:
    doc: BAM files with marked duplicate reads.
    type: File[]
    outputSource: index_dups_marked_bams/indexed_file
  output_pbc_files:
    doc: PCR Bottleneck Coeficient files.
    type: File[]
    outputSource: execute_pcr_bottleneck_coef/pbc_file
  output_percent_mitochondrial_reads:
    doc: Percentage of mitochondrial reads.
    type: File[]
    outputSource: percent_mitochondrial_reads/percent_map
  output_percentage_uniq_reads:
    doc: Percentage of uniq reads from preseq c_curve output
    type: File[]
    outputSource: percent_uniq_reads/output
  output_picard_mark_duplicates_files:
    doc: Picard MarkDuplicates metrics files.
    type: File[]
    outputSource: mark_duplicates/output_metrics_file
  output_preseq_c_curve_files:
    doc: Preseq c_curve output files.
    type: File[]
    outputSource: preseq-c-curve/output_file
  output_read_count_mapped:
    doc: Read counts of the mapped BAM files
    type: File[]
    outputSource: mapped_reads_count/output
  output_read_count_mapped_filtered:
    doc: Read counts of the mapped and filtered BAM files
    type: File[]
    outputSource: mapped_filtered_reads_count/output_read_count

steps:
  bam_idxstats:
    in:
      bam: index_bams/indexed_file
    scatter: bam
    run: map.cwl.steps/bam_idxstats.cwl
    out:
    - idxstats_file
  bowtie-se:
    in:
      X:
        valueFrom: ${return 2000}
      genome_ref_first_index_file: genome_ref_first_index_file
      input_fastq_file: input_fastq_files
      nthreads: nthreads
      output_filename: extract_basename_2/output_path
      v:
        valueFrom: ${return 2}
    scatter:
    - input_fastq_file
    - output_filename
    scatterMethod: dotproduct
    run: map.cwl.steps/bowtie-se.cwl
    out:
    - output_aligned_file
    - output_bowtie_log
  execute_pcr_bottleneck_coef:
    in:
      genome_sizes: genome_sizes_file
      input_bam_files: filtered2sorted/sorted_file
      input_output_filenames: extract_basename_2/output_path
    run: map.cwl.steps/execute_pcr_bottleneck_coef.cwl
    out:
    - pbc_file
  extract_basename_1:
    in:
      input_file: input_fastq_files
    scatter: input_file
    run: map.cwl.steps/extract_basename_1.cwl
    out:
    - output_basename
  extract_basename_2:
    in:
      file_path: extract_basename_1/output_basename
    scatter: file_path
    run: map.cwl.steps/extract_basename_2.cwl
    out:
    - output_path
  filter-unmapped:
    in:
      input_file: sort_bams/sorted_file
      output_filename: extract_basename_2/output_path
    scatter:
    - input_file
    - output_filename
    scatterMethod: dotproduct
    run: map.cwl.steps/filter-unmapped.cwl
    out:
    - filtered_file
  filtered2sorted:
    in:
      input_file: filter-unmapped/filtered_file
      nthreads: nthreads
    scatter:
    - input_file
    run: map.cwl.steps/filtered2sorted.cwl
    out:
    - sorted_file
  index_bams:
    in:
      input_file: sort_bams/sorted_file
    scatter: input_file
    run: map.cwl.steps/index_bams.cwl
    out:
    - indexed_file
  index_dedup_bams:
    in:
      input_file: sort_dedup_bams/sorted_file
    scatter:
    - input_file
    run: map.cwl.steps/index_dedup_bams.cwl
    out:
    - indexed_file
  index_dups_marked_bams:
    in:
      input_file: sort_dups_marked_bams/sorted_file
    scatter:
    - input_file
    run: map.cwl.steps/index_dups_marked_bams.cwl
    out:
    - indexed_file
  index_filtered_bam:
    in:
      input_file: filtered2sorted/sorted_file
    scatter: input_file
    run: map.cwl.steps/index_filtered_bam.cwl
    out:
    - indexed_file
  mapped_filtered_reads_count:
    in:
      input_bam_file: sort_dedup_bams/sorted_file
      output_suffix:
        valueFrom: .mapped_and_filtered.read_count.txt
    scatter: input_bam_file
    run: map.cwl.steps/mapped_filtered_reads_count.cwl
    out:
    - output_read_count
  mapped_reads_count:
    in:
      bowtie_log: bowtie-se/output_bowtie_log
    scatter: bowtie_log
    run: map.cwl.steps/mapped_reads_count.cwl
    out:
    - output
  mark_duplicates:
    in:
      input_file: index_filtered_bam/indexed_file
      java_opts: picard_java_opts
      output_filename: extract_basename_2/output_path
      output_suffix:
        valueFrom: bam
      picard_jar_path: picard_jar_path
    scatter:
    - input_file
    - output_filename
    scatterMethod: dotproduct
    run: map.cwl.steps/mark_duplicates.cwl
    out:
    - output_metrics_file
    - output_dedup_bam_file
  percent_mitochondrial_reads:
    in:
      chrom:
        valueFrom: chrM
      idxstats: bam_idxstats/idxstats_file
      output_filename:
        valueFrom: |-
          ${return inputs.idxstats.basename.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '').replace(/\.[^/.]+$/, '').replace(/\.[^/.]+$/, '.mitochondrial_percentage.txt')}
    scatter: idxstats
    run: map.cwl.steps/percent_mitochondrial_reads.cwl
    out:
    - percent_map
  percent_uniq_reads:
    in:
      preseq_c_curve_outfile: preseq-c-curve/output_file
    scatter: preseq_c_curve_outfile
    run: map.cwl.steps/percent_uniq_reads.cwl
    out:
    - output
  preseq-c-curve:
    in:
      input_sorted_file: filtered2sorted/sorted_file
      output_file_basename: extract_basename_2/output_path
    scatter:
    - input_sorted_file
    - output_file_basename
    scatterMethod: dotproduct
    run: map.cwl.steps/preseq-c-curve.cwl
    out:
    - output_file
  remove_duplicates:
    in:
      F:
        valueFrom: ${return 1024}
      b:
        valueFrom: ${return true}
      input_file: index_dups_marked_bams/indexed_file
      outfile_name:
        valueFrom: ${return inputs.input_file.basename.replace('dups_marked', 'dedup')}
      suffix:
        valueFrom: .dedup.bam
    scatter:
    - input_file
    run: map.cwl.steps/remove_duplicates.cwl
    out:
    - outfile
  sam2bam:
    in:
      input_file: bowtie-se/output_aligned_file
      nthreads: nthreads
    scatter: input_file
    run: map.cwl.steps/sam2bam.cwl
    out:
    - bam_file
  sort_bams:
    in:
      input_file: sam2bam/bam_file
      nthreads: nthreads
    scatter: input_file
    run: map.cwl.steps/sort_bams.cwl
    out:
    - sorted_file
  sort_dedup_bams:
    in:
      input_file: remove_duplicates/outfile
      nthreads: nthreads
    scatter:
    - input_file
    run: map.cwl.steps/sort_dedup_bams.cwl
    out:
    - sorted_file
  sort_dups_marked_bams:
    in:
      input_file: mark_duplicates/output_dedup_bam_file
      nthreads: nthreads
      suffix:
        valueFrom: .dups_marked.bam
    scatter:
    - input_file
    run: map.cwl.steps/sort_dups_marked_bams.cwl
    out:
    - sorted_file
