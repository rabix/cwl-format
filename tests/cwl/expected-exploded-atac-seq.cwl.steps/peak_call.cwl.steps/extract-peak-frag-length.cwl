#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Extracts best fragment length from SPP output text file

inputs:
  input_spp_txt_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_best_frag_length:
    type: float
    outputBinding:
      glob: best_frag_length
      outputEval: $(Number(self[0].contents.replace('\n', '')))
      loadContents: true
stdout: best_frag_length

baseCommand: extract-best-frag-length.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
