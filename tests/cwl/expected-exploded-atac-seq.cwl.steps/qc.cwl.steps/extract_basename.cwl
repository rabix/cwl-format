#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Extracts the base name of a file

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_basename:
    type: string
    outputBinding:
      outputEval: |-
        $(inputs.input_file.path.substr(inputs.input_file.path.lastIndexOf('/') + 1, inputs.input_file.path.lastIndexOf('.') - (inputs.input_file.path.lastIndexOf('/') + 1)))

baseCommand: echo

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
