#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Filter BAM file to only include reads overlapping with a BED file

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_bam_file:
    doc: Aligned BAM file to filter
    type: File
    inputBinding:
      position: 3
  input_bedfile:
    doc: Bedfile used to only include reads overlapping this BED FILE
    type: File
    inputBinding:
      prefix: -L
      position: 2

outputs:
  filtered_file:
    doc: Filtered aligned BAM file by BED coordinates file
    type: File
    outputBinding:
      glob: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '') + '.in_peaks.bam')
stdout: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '') + '.in_peaks.bam')

baseCommand:
- samtools
- view
- -b
- -h

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
