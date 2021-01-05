#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |-
  Usage: c_curve [OPTIONS] <sorted-bed-file>

  Options:
    -o, -output   yield output file (default: stdout) 
    -s, -step     step size in extrapolations (default: 1e+06) 
    -v, -verbose  print more information 
    -P, -pe       input is paired end read file 
    -H, -hist     input is a text file containing the observed histogram 
    -V, -vals     input is a text file containing only the observed counts 
    -B, -bam      input is in BAM format 
    -l, -seg_len  maximum segment length when merging paired end bam reads 
                  (default: 5000) 

  Help options:
    -?, -help     print this help message 
        -about    print about message

requirements:
  InlineJavascriptRequirement: {}

inputs:
  B:
    doc: "-bam      input is in BAM format \n"
    type: boolean
    default: true
    inputBinding:
      prefix: -B
      position: 1
  H:
    doc: "-hist     input is a text file containing the observed histogram \n"
    type: File?
    inputBinding:
      prefix: -H
      position: 1
  V:
    doc: "-vals     input is a text file containing only the observed counts \n"
    type: File?
    inputBinding:
      prefix: -V
      position: 1
  input_sorted_file:
    doc: Sorted bed or BAM file
    type: File
    inputBinding:
      position: 2
  l:
    doc: |
      -seg_len  maximum segment length when merging paired end bam reads 
      (default: 5000)
      Help options:
      -?, -help     print this help message
      -about    print about message
    type: int?
    inputBinding:
      prefix: -l
      position: 1
  output_file_basename:
    type: string
  pe:
    doc: "-pe       input is paired end read file \n"
    type: boolean?
    inputBinding:
      prefix: -P
      position: 1
  s:
    doc: "-step     step size in extrapolations (default: 1e+06) \n"
    type: float?
    inputBinding:
      prefix: -s
      position: 1
  v:
    doc: "-verbose  print more information \n"
    type: boolean
    default: false
    inputBinding:
      prefix: -v
      position: 1

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output_file_basename + '.preseq_c_curve.txt')
stdout: $(inputs.output_file_basename + '.preseq_c_curve.txt')

baseCommand:
- preseq
- c_curve

hints:
  DockerRequirement:
    dockerPull: reddylab/preseq:2.0
