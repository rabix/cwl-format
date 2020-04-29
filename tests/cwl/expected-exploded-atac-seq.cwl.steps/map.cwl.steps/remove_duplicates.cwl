class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}

inputs:
  F:
    doc: only include reads with none of the bits set in INT set in FLAG [0]
    type: int?
    inputBinding:
      prefix: -F
      position: 1
  L:
    doc: FILE  only include reads overlapping this BED FILE [null]
    type: File?
    inputBinding:
      prefix: -L
      position: 1
  S:
    doc: Input format autodetected
    type: boolean
    default: true
    inputBinding:
      prefix: -S
      position: 1
  b:
    doc: output BAM
    type: boolean?
    inputBinding:
      prefix: -b
      position: 1
  f:
    doc: only include reads with all bits set in INT set in FLAG [0]
    type: int?
    inputBinding:
      prefix: -f
      position: 1
  header:
    doc: Include header in output
    type: boolean?
    inputBinding:
      prefix: -h
      position: 1
  input_file:
    doc: File to be converted to BAM with samtools
    type: File
    inputBinding:
      position: 2
  nthreads:
    doc: Number of threads used
    type: int
    default: 1
    inputBinding:
      prefix: -@
      position: 1
  outfile_name:
    doc: |-
      Output file name. If not specified, the basename of the input file with the suffix specified in the suffix argument will be used.
    type: string?
  q:
    doc: only include reads with mapping quality >= INT [0]
    type: int?
    inputBinding:
      prefix: -q
      position: 1
  suffix:
    doc: suffix of the transformed SAM/BAM file (including extension, e.g. .filtered.sam)
    type: string?
  u:
    doc: uncompressed BAM output (implies -b)
    type: boolean
    default: true
    inputBinding:
      prefix: -u
      position: 1

outputs:
  outfile:
    doc: Aligned file in SAM or BAM format
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outfile_name) return inputs.outfile_name;
          var suffix = inputs.b ? '.bam' : '.sam';
          suffix = inputs.suffix || suffix;
          return inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + suffix
        }
stdout: |
  ${
     if (inputs.outfile_name) return inputs.outfile_name;
     var suffix = inputs.b ? '.bam' : '.sam';
     suffix = inputs.suffix || suffix;
     return inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + suffix
  }

baseCommand:
- samtools
- view

hints:
  DockerRequirement:
    dockerPull: dukegcb/samtools:1.3
