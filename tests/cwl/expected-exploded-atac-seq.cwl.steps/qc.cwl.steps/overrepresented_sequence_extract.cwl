class: CommandLineTool
cwlVersion: v1.0

inputs:
  default_adapters_file:
    doc: Adapters file in fasta format
    type: File
    inputBinding:
      position: 2
  input_basename:
    doc: Name of the sample - used as a base name for generating output files
    type: string
  input_fastqc_data:
    doc: fastqc_data.txt file from a fastqc report
    type: File
    inputBinding:
      position: 1

outputs:
  output_custom_adapters:
    type: File
    outputBinding:
      glob: $(inputs.input_basename + '.custom_adapters.fasta')

baseCommand: overrepresented_sequence_extract.py
arguments:
- position: 3
  valueFrom: $(inputs.input_basename + '.custom_adapters.fasta')

hints:
  DockerRequirement:
    dockerPull: reddylab/overrepresented_sequence_extract:1.0
