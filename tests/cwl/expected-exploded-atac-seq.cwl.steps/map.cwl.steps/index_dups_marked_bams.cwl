class: CommandLineTool
cwlVersion: v1.0

requirements:
  InitialWorkDirRequirement:
    listing:
    - $(inputs.input_file)
  InlineJavascriptRequirement: {}

inputs:
  input_file:
    doc: Aligned file to be sorted with samtools
    type: File
    inputBinding:
      position: 1

outputs:
  indexed_file:
    doc: Indexed BAM file
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.input_file.basename)

baseCommand:
- samtools
- index
arguments:
- position: 2
  valueFrom: $(inputs.input_file.basename + '.bai')

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
