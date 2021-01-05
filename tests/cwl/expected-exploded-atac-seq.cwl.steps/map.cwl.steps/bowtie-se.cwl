#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  X:
    doc: 'maximum insert size for paired-end alignment (default: 250)'
    type: int?
    inputBinding:
      prefix: -X
      position: 4
  best:
    doc: Hits guaranteed best stratum; ties broken by quality
    type: boolean
    default: true
    inputBinding:
      prefix: --best
      position: 5
  chunkmbs:
    doc: |-
      The number of megabytes of memory a given thread is given to store path descriptors in --best mode. (Default: 256)
    type: int?
    inputBinding:
      prefix: --chunkmbs
      position: 5
  genome_ref_first_index_file:
    doc: |-
      First file (extension .1.ebwt) of the Bowtie index files generated for the reference genome (see http://bowtie-bio.sourceforge.net/tutorial.shtml#newi)
    type: File
    secondaryFiles:
    - ^^.2.ebwt
    - ^^.3.ebwt
    - ^^.4.ebwt
    - ^^.rev.1.ebwt
    - ^^.rev.2.ebwt
    inputBinding:
      position: 9
      valueFrom: $(self.path.split('.').splice(0,self.path.split('.').length-2).join("."))
  input_fastq_file:
    doc: Query input FASTQ file.
    type: File
    inputBinding:
      position: 10
  m:
    doc: 'Suppress all alignments if > <int> exist (def: 1)'
    type: int
    default: 1
    inputBinding:
      prefix: -m
      position: 7
  nthreads:
    doc: '<int> number of alignment threads to launch (default: 1)'
    type: int
    default: 1
    inputBinding:
      prefix: --threads
      position: 8
  output_filename:
    type: string
  sam:
    doc: 'Write hits in SAM format (default: BAM)'
    type: boolean
    default: true
    inputBinding:
      prefix: --sam
      position: 2
  seedlen:
    doc: 'seed length for -n (default: 28)'
    type: int?
    inputBinding:
      prefix: --seedlen
      position: 1
  seedmms:
    doc: 'max mismatches in seed (between [0, 3], default: 2)'
    type: int?
    inputBinding:
      prefix: --seedmms
      position: 1
  strata:
    doc: Hits in sub-optimal strata aren't reported (requires --best)
    type: boolean
    default: true
    inputBinding:
      prefix: --strata
      position: 6
  t:
    doc: Print wall-clock time taken by search phases
    type: boolean
    default: true
    inputBinding:
      prefix: -t
      position: 1
  trim3:
    doc: trim <int> bases from 3' (right) end of reads
    type: int?
    inputBinding:
      prefix: --trim3
      position: 1
  trim5:
    doc: trim <int> bases from 5' (left) end of reads
    type: int?
    inputBinding:
      prefix: --trim5
      position: 1
  v:
    doc: Report end-to-end hits w/ <=v mismatches; ignore qualities
    type: int?
    inputBinding:
      prefix: -v
      position: 3

outputs:
  output_aligned_file:
    doc: Aligned bowtie file in [SAM|BAM] format.
    type: File
    outputBinding:
      glob: $(inputs.output_filename + '.sam')
  output_bowtie_log:
    type: File
    outputBinding:
      glob: $(inputs.output_filename + '.bowtie.log')
stderr: $(inputs.output_filename + '.bowtie.log')

baseCommand: bowtie
arguments:
- position: 11
  valueFrom: $(inputs.output_filename + '.sam')

hints:
  DockerRequirement:
    dockerPull: dukegcb/bowtie
