#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Counts lines in a file and returns a suffixed file with that number

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_file:
    type: File
  output_suffix:
    type: string
    default: .count

outputs:
  output_counts:
    type: File
    outputBinding:
      glob: $(inputs.input_file.path.replace(/^.*[\\\/]/, '') + inputs.output_suffix)
stdout: $(inputs.input_file.path.replace(/^.*[\\\/]/, '') + inputs.output_suffix)

baseCommand:
- wc
- -l
stdin: $(inputs.input_file.path)
