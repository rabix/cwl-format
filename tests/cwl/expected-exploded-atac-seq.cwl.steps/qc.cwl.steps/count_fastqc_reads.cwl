#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Extracts read count from fastqc_data.txt

inputs:
  input_basename:
    type: string
  input_fastqc_data:
    type: File
    inputBinding:
      position: 1

outputs:
  output_fastqc_read_count:
    type: File
    outputBinding:
      glob: $(inputs.input_basename + '.fastqc-read_count.txt')
stdout: $(inputs.input_basename + '.fastqc-read_count.txt')

baseCommand: count-fastqc_data-reads.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
