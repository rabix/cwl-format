#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Get number of processed reads from Bowtie log.

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  preseq_c_curve_outfile:
    type: File
    inputBinding: {}

outputs:
  output:
    type: File
    outputBinding:
      glob: |-
        $(inputs.preseq_c_curve_outfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + '.percentage_unique_reads.txt')
stdout: |-
  $(inputs.preseq_c_curve_outfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + '.percentage_unique_reads.txt')

baseCommand: percent-uniq-reads-from-preseq.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
