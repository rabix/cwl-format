
# Top comment is preserved

class: CommandLineTool
cwlVersion: v1.0
requirements:
  InlineJavascriptRequirement: {}
inputs:
  in1:
    type: string
    inputBinding:
      position: 1
      valueFrom: A_$(inputs.in1)_B_${return inputs.in1}_C_$(inputs.in1)
outputs:
  out1:
    type: string
    outputBinding:
      glob: out.txt
      outputEval: $(self[0].contents)_D_$(runtime.cores)
      loadContents: true
baseCommand: echo
arguments:
- valueFrom: $(runtime)
stdout: out.txt
