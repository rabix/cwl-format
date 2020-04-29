class: Workflow
cwlVersion: v1.0
doc: ATAC-seq 04 quantification - SE

requirements:
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement

inputs:
  as_narrowPeak_file:
    doc: Definition narrowPeak file in AutoSql format (used in bedToBigBed)
    type: File
  genome_effective_size:
    doc: |-
      Effective genome size used by MACS2. It can be numeric or a shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
    type: string
    default: hs
  input_bam_files:
    type: File[]
  input_genome_sizes:
    doc: Two column tab-delimited file with chromosome size information
    type: File
  nthreads:
    type: int
    default: 1

outputs:
  output_extended_peak_file:
    doc: peakshift/phantomPeak extended fragment results file
    type: File[]
    outputSource: peak-calling/output_ext_frag_bdg_file
  output_filtered_read_count_file:
    doc: Filtered read count reported by MACS2
    type: File[]
    outputSource: count-reads-filtered/read_count_file
  output_peak_bigbed_file:
    doc: Peaks in bigBed format
    type: File[]
    outputSource: peaks-bed-to-bigbed/bigbed
  output_peak_count_within_replicate:
    doc: Peak counts within replicate
    type: File[]
    outputSource: count-peaks/output_counts
  output_peak_file:
    doc: peakshift/phantomPeak results file
    type: File[]
    outputSource: peak-calling/output_peak_file
  output_peak_summits_file:
    doc: File containing peak summits
    type: File[]
    outputSource: peak-calling/output_peak_summits_file
  output_peak_xls_file:
    doc: Peak calling report file (*_peaks.xls file produced by MACS2)
    type: File[]
    outputSource: peak-calling/output_peak_xls_file
  output_read_in_peak_count_within_replicate:
    doc: Reads peak counts within replicate
    type: File[]
    outputSource: extract-count-reads-in-peaks/output_read_count
  output_spp_cross_corr_plot:
    doc: peakshift/phantomPeak results file
    type: File[]
    outputSource: spp/output_spp_cross_corr_plot
  output_spp_x_cross_corr:
    doc: peakshift/phantomPeak results file
    type: File[]
    outputSource: spp/output_spp_cross_corr

steps:
  count-peaks:
    in:
      input_file: peak-calling/output_peak_file
      output_suffix:
        valueFrom: .peak_count.within_replicate.txt
    scatter: input_file
    run: peak_call.cwl.steps/count-peaks.cwl
    out:
    - output_counts
  count-reads-filtered:
    in:
      peak_xls_file: peak-calling/output_peak_xls_file
    scatter: peak_xls_file
    run: peak_call.cwl.steps/count-reads-filtered.cwl
    out:
    - read_count_file
  extract-count-reads-in-peaks:
    in:
      input_bam_file: filter-reads-in-peaks/filtered_file
      output_suffix:
        valueFrom: .read_count.within_replicate.txt
    scatter: input_bam_file
    run: peak_call.cwl.steps/extract-count-reads-in-peaks.cwl
    out:
    - output_read_count
  extract-peak-frag-length:
    in:
      input_spp_txt_file: spp/output_spp_cross_corr
    scatter: input_spp_txt_file
    run: peak_call.cwl.steps/extract-peak-frag-length.cwl
    out:
    - output_best_frag_length
  filter-reads-in-peaks:
    in:
      input_bam_file: input_bam_files
      input_bedfile: peak-calling/output_peak_file
    scatter:
    - input_bam_file
    - input_bedfile
    scatterMethod: dotproduct
    run: peak_call.cwl.steps/filter-reads-in-peaks.cwl
    out:
    - filtered_file
  peak-calling:
    in:
      format:
        valueFrom: BAM
      bdg:
        valueFrom: ${return true}
      extsize:
        valueFrom: ${return 200}
      g: genome_effective_size
      nomodel:
        valueFrom: ${return true}
      q:
        valueFrom: ${return 0.1}
      shift:
        valueFrom: ${return -100}
      treatment:
        valueFrom: $([self])
        source: input_bam_files
    scatter:
    - treatment
    scatterMethod: dotproduct
    run: peak_call.cwl.steps/peak-calling.cwl
    out:
    - output_peak_file
    - output_peak_summits_file
    - output_ext_frag_bdg_file
    - output_peak_xls_file
  peaks-bed-to-bigbed:
    in:
      type:
        valueFrom: bed6+4
      as: as_narrowPeak_file
      bed: trunk-peak-score/trunked_scores_peaks
      genome_sizes: input_genome_sizes
    scatter: bed
    run: peak_call.cwl.steps/peaks-bed-to-bigbed.cwl
    out:
    - bigbed
  spp:
    in:
      input_bam: input_bam_files
      nthreads: nthreads
      savp:
        valueFrom: ${return true}
    scatter:
    - input_bam
    scatterMethod: dotproduct
    run: peak_call.cwl.steps/spp.cwl
    out:
    - output_spp_cross_corr
    - output_spp_cross_corr_plot
  trunk-peak-score:
    in:
      peaks: peak-calling/output_peak_file
    scatter: peaks
    run: peak_call.cwl.steps/trunk-peak-score.cwl
    out:
    - trunked_scores_peaks
