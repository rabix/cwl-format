#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

inputs:
  S:
    doc: Input format autodetected
    type: boolean
    default: true
    inputBinding:
      prefix: -S
      position: 1
  input_file:
    doc: File to be converted to BAM with samtools
    type: File
    inputBinding:
      position: 2
  nthreads:
    doc: Number of threads used
    type: int
    default: 1
    inputBinding:
      prefix: -@
      position: 1

outputs:
  bam_file:
    doc: Aligned file in BAM format
    type: File
    outputBinding:
      glob: |-
        $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.bam')
stdout: |-
  $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.bam')

baseCommand:
- samtools
- view
- -b

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
