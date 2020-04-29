class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_file:
    doc: Aligned file to be sorted with samtools
    type: File
    inputBinding:
      position: 1000
  n:
    doc: Sort by read name
    type: boolean
    default: false
    inputBinding:
      prefix: -n
      position: 1
  nthreads:
    doc: Number of threads used in sorting
    type: int
    default: 1
    inputBinding:
      prefix: -@
      position: 1
  suffix:
    doc: suffix of the transformed SAM/BAM file (including extension, e.g. .filtered.sam)
    type: string
    default: .sorted.bam

outputs:
  sorted_file:
    doc: Sorted aligned file
    type: File
    outputBinding:
      glob: |-
        $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + inputs.suffix)
stdout: |-
  $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + inputs.suffix)

baseCommand:
- samtools
- sort

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
