#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_file:
    doc: Aligned file to be sorted with samtools
    type: File
    inputBinding:
      position: 1000
  nthreads:
    doc: Number of threads used in sorting
    type: int
    default: 1
    inputBinding:
      prefix: -@
      position: 1
  output_filename:
    doc: Basename for the output file
    type: string

outputs:
  filtered_file:
    doc: Filter unmapped reads in aligned file
    type: File
    outputBinding:
      glob: $(inputs.output_filename + '.accepted_hits.bam')
stdout: $(inputs.output_filename + '.accepted_hits.bam')

baseCommand:
- samtools
- view
- -F
- '4'
- -b
- -h

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
