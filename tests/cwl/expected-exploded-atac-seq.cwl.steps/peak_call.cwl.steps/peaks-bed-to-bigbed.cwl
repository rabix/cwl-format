#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: |
  "bedToBigBed v. 2.7 - Convert bed file to bigBed. (BigBed version: 4)
  usage:
     bedToBigBed in.bed chrom.sizes out.bb
  Where in.bed is in one of the ascii bed formats, but not including track lines
  and chrom.sizes is two column: <chromosome name> <size in bases>
  and out.bb is the output indexed big bed file.
  Use the script: fetchChromSizes to obtain the actual chrom.sizes information
  from UCSC, please do not make up a chrom sizes from your own information.
  The in.bed file must be sorted by chromosome,start,
    to sort a bed file, use the unix sort command:
       sort -k1,1 -k2,2n unsorted.bed > sorted.bed"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  type:
    doc: |
      -type=bedN[+[P]] :
                      N is between 3 and 15,
                      optional (+) if extra "bedPlus" fields,
                      optional P specifies the number of extra fields. Not required, but preferred.
                      Examples: -type=bed6 or -type=bed6+ or -type=bed6+3
                      (see http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    type: string?
    inputBinding:
      prefix: -type=
      position: 1
      separate: false
  as:
    doc: |
      -as=fields.as - If you have non-standard "bedPlus" fields, it's great to put a definition
                       of each field in a row in AutoSql format here.1)
    type: File?
    inputBinding:
      prefix: -as=
      position: 1
      separate: false
  bed:
    doc: Input bed file
    type: File
    inputBinding:
      position: 2
  blockSize:
    doc: "-blockSize=N - Number of items to bundle in r-tree.  Default 256\n"
    type: int?
    inputBinding:
      prefix: -blockSize=
      position: 1
      separate: false
  extraIndex:
    doc: |
      -extraIndex=fieldList - If set, make an index on each field in a comma separated list
         extraIndex=name and extraIndex=name,id are commonly used.
    type:
    - 'null'
    - type: array
      items: string
    inputBinding:
      prefix: -extraIndex=
      position: 1
      itemSeparator: ','
  genome_sizes:
    doc: "genome_sizes is two column: <chromosome name> <size in bases>.\n"
    type: File
    inputBinding:
      position: 3
  itemsPerSlot:
    doc: "-itemsPerSlot=N - Number of data points bundled at lowest level. Default\
      \ 512\n"
    type: int?
    inputBinding:
      prefix: -itemsPerSlot=
      position: 1
      separate: false
  output_suffix:
    type: string
    default: .bb
  tab:
    doc: |
      -tab - If set, expect fields to be tab separated, normally expects white space separator.
    type: boolean?
    inputBinding:
      position: 1
  unc:
    doc: "-unc - If set, do not use compression.\n"
    type: boolean?
    inputBinding:
      position: 1

outputs:
  bigbed:
    type: File
    outputBinding:
      glob: $(inputs.bed.path.replace(/^.*[\\\/]/, '')+ inputs.output_suffix)

baseCommand: bedToBigBed
arguments:
- position: 4
  valueFrom: $(inputs.bed.path.replace(/^.*[\\\/]/, '') + inputs.output_suffix)

hints:
  DockerRequirement:
    dockerPull: dleehr/docker-hubutils
