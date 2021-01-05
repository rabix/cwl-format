#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  bedSort - Sort a .bed file by chrom,chromStart
  usage:
     bedSort in.bed out.bed
  in.bed and out.bed may be the same.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  bed_file:
    doc: Bed or bedGraph file to be sorted
    type: File
    inputBinding:
      position: 1

outputs:
  bed_file_sorted:
    type: File
    outputBinding:
      glob: $(inputs.bed_file.path.replace(/^.*[\\\/]/, '') + "_sorted")

baseCommand: bedSort
arguments:
- position: 2
  valueFrom: $(inputs.bed_file.path.replace(/^.*[\\\/]/, '') + "_sorted")

hints:
  DockerRequirement:
    dockerPull: dleehr/docker-hubutils
