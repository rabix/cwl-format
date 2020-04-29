class: CommandLineTool
cwlVersion: v1.0
doc: Count number of dedup-ed reads used in peak calling

requirements:
  InlineJavascriptRequirement: {}

inputs:
  peak_xls_file:
    type: File
    inputBinding:
      position: 1

outputs:
  read_count_file:
    type: File
    outputBinding:
      glob: |-
        $(inputs.peak_xls_file.path.replace(/^.*[\\\/]/, '').replace(/\_peaks\.xls$/, '_read_count.txt'))
stdout: |-
  $(inputs.peak_xls_file.path.replace(/^.*[\\\/]/, '').replace(/\_peaks\.xls$/, '_read_count.txt'))

baseCommand: count-filtered-reads-macs2.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
