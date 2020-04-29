class: CommandLineTool
cwlVersion: v1.0
doc: Compares 2 files

inputs:
  brief:
    type: boolean
    default: true
    inputBinding:
      prefix: --brief
      position: 3
  file1:
    type: File
    inputBinding:
      position: 1
  file2:
    type: File
    inputBinding:
      position: 2

outputs:
  result:
    type: File
    outputBinding:
      glob: stdout.txt
stdout: stdout.txt

baseCommand: diff

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
