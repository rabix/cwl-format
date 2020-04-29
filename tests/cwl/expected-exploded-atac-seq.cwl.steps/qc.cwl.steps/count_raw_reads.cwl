class: CommandLineTool
cwlVersion: v1.0
doc: Counts reads in a fastq file

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_basename:
    type: string
  input_fastq_file:
    type: File
    inputBinding:
      position: 1

outputs:
  output_read_count:
    type: File
    outputBinding:
      glob: $(inputs.input_basename + '.read_count.txt')
stdout: $(inputs.input_basename + '.read_count.txt')

baseCommand: count-fastq-reads.sh

hints:
  DockerRequirement:
    dockerPull: reddylab/workflow-utils:ggr
