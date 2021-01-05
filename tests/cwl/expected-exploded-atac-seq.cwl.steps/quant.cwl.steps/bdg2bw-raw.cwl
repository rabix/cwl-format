#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: 'Tool:   bedGraphToBigWig v 4 - Convert a bedGraph file to bigWig format.'

requirements:
  InlineJavascriptRequirement: {}

inputs:
  bed_graph:
    doc: "\tbed_graph is a four column file in the format: <chrom> <start> <end> <value>\n"
    type: File
    inputBinding:
      position: 1
  genome_sizes:
    doc: "\tgenome_sizes is two column: <chromosome name> <size in bases>.\n"
    type: File
    inputBinding:
      position: 2
  output_suffix:
    type: string
    default: .bw

outputs:
  output_bigwig:
    type: File
    outputBinding:
      glob: |-
        $(inputs.bed_graph.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

baseCommand: bedGraphToBigWig
arguments:
- position: 3
  valueFrom: |-
    $(inputs.bed_graph.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

hints:
  DockerRequirement:
    dockerPull: dukegcb/bedgraphtobigwig
