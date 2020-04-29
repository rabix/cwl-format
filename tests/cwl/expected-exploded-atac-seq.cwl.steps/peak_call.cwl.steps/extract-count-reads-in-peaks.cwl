class: CommandLineTool
cwlVersion: v1.0
doc: Extract mapped reads from BAM file using Samtools flagstat command

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  input_bam_file:
    doc: Aligned BAM file to filter
    type: File
    inputBinding:
      position: 1
  output_suffix:
    type: string

outputs:
  output_read_count:
    doc: Samtools Flagstat report file
    type: File
    outputBinding:
      glob: |-
        $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)
stdout: |-
  $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

baseCommand:
- samtools
- flagstat
arguments:
- position: 10000
  valueFrom: " | head -n1 | cut -f 1 -d ' '"
  shellQuote: false

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools
