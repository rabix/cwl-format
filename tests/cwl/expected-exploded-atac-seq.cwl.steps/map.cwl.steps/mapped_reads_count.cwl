#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Get number of processed reads from Bowtie log.

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  bowtie_log:
    type: File
    inputBinding: {}

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.bowtie_log.path.replace(/^.*[\\\/]/, '') + '.read_count.mapped')
stdout: $(inputs.bowtie_log.path.replace(/^.*[\\\/]/, '') + '.read_count.mapped')

baseCommand: read-count-from-bowtie-log.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
