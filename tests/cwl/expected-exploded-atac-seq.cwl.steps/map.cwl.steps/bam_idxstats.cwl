class: CommandLineTool
cwlVersion: v1.0

requirements:
  InitialWorkDirRequirement:
    listing:
    - $(inputs.bam)
  InlineJavascriptRequirement: {}

inputs:
  bam:
    doc: Bam file (it should be indexed)
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 1

outputs:
  idxstats_file:
    doc: |
      Idxstats output file. TAB-delimited with each line consisting of reference
      sequence name, sequence length, # mapped reads and # unmapped reads
    type: File
    outputBinding:
      glob: $(inputs.bam.basename + ".idxstats")
stdout: $(inputs.bam.basename + ".idxstats")

baseCommand:
- samtools
- idxstats

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
