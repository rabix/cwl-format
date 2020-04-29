class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

inputs:
  format:
    type: string
    default: fastq
    inputBinding:
      prefix: --format
      position: 3
  input_fastq_file:
    type: File
    inputBinding:
      position: 4
  noextract:
    type: boolean
    default: true
    inputBinding:
      prefix: --noextract
      position: 2
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --threads
      position: 5

outputs:
  output_qc_report_file:
    type: File
    outputBinding:
      glob: |-
        $(inputs.input_fastq_file.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, '') + "_fastqc.zip")

baseCommand: fastqc
arguments:
- prefix: --dir
  position: 5
  valueFrom: $(runtime.tmpdir)
- prefix: -o
  position: 5
  valueFrom: $(runtime.outdir)

hints:
  DockerRequirement:
    dockerPull: dukegcb/fastqc
