class: CommandLineTool
cwlVersion: v1.0
doc: |-
  Unzips a zipped fastqc report and returns the fastqc_data.txt file. Unzips the file to pipe and uses redirection

inputs:
  extract_pattern:
    type: string
    default: '*/fastqc_data.txt'
    inputBinding:
      position: 3
  input_basename:
    type: string
  input_qc_report_file:
    type: File
    inputBinding:
      position: 2
  pipe:
    type: boolean
    default: true
    inputBinding:
      prefix: -p
      position: 1

outputs:
  output_fastqc_data_file:
    type: File
    outputBinding:
      glob: $(inputs.input_basename + '.fastqc_data.txt')
stdout: $(inputs.input_basename + '.fastqc_data.txt')

baseCommand: unzip

hints:
  DockerRequirement:
    dockerPull: dukegcb/fastqc
