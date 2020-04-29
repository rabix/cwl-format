class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  barcode_tag:
    doc: |-
      If true do not write duplicates to the output file instead of writing them with appropriate flags set.  (Default false).
    type: string?
    inputBinding:
      prefix: BARCODE_TAG=
      position: 5
      separate: false
  input_file:
    doc: One or more input SAM or BAM files to analyze. Must be coordinate sorted.
    type: File
    inputBinding:
      position: 4
      valueFrom: $('INPUT=' + self.path)
      shellQuote: false
  java_opts:
    doc: |-
      JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
    type: string?
    inputBinding:
      position: 1
      shellQuote: false
  metrics_suffix:
    doc: 'Suffix used to create the metrics output file (Default: dedup_metrics.txt)'
    type: string
    default: dedup_metrics.txt
  output_filename:
    doc: Output filename used as basename
    type: string
  output_suffix:
    doc: 'Suffix used to identify the output file (Default: dedup.bam)'
    type: string
    default: dedup.bam
  picard_jar_path:
    doc: Path to the picard.jar file
    type: string
    inputBinding:
      prefix: -jar
      position: 2
  remove_duplicates:
    doc: |-
      If true do not write duplicates to the output file instead of writing them with appropriate flags set.  (Default false).
    type: boolean
    default: false
    inputBinding:
      position: 5
      valueFrom: $('REMOVE_DUPLICATES=' + self)

outputs:
  output_dedup_bam_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename  + '.' + inputs.output_suffix)
  output_metrics_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename + '.' + inputs.metrics_suffix)

baseCommand:
- java
arguments:
- position: 3
  valueFrom: MarkDuplicates
- position: 5
  valueFrom: $('OUTPUT=' + inputs.output_filename + '.' + inputs.output_suffix)
  shellQuote: false
- position: 5
  valueFrom: $('METRICS_FILE='+inputs.output_filename + '.' + inputs.metrics_suffix)
  shellQuote: false
- position: 5
  valueFrom: $('TMP_DIR='+runtime.tmpdir)
  shellQuote: false

hints:
  DockerRequirement:
    dockerPull: dukegcb/picard
