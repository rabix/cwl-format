class: CommandLineTool
cwlVersion: v1.0

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  control_bam:
    doc: |-
      <Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
    type: File?
    inputBinding:
      prefix: -i=
      separate: false
  fdr:
    doc: -fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
    type: float?
    inputBinding:
      prefix: -fdr=
      separate: false
  filtchr:
    doc: |-
      -filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name
    type: string?
    inputBinding:
      prefix: -filtchr=
      separate: false
  input_bam:
    doc: |-
      <ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped)(FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
    type: File
    inputBinding:
      prefix: -c=
      separate: false
  npeak:
    doc: -npeak=<numPeaks>, threshold on number of peaks to call
    type: int?
    inputBinding:
      prefix: -npeak=
      separate: false
  nthreads:
    doc: -p=<nodes> , number of parallel processing nodes, default=0
    type: int?
    inputBinding:
      prefix: -p=
      separate: false
  rf:
    doc: 'overwrite (force remove) output files in case they exist. Default: true'
    type: boolean
    default: true
    inputBinding:
      prefix: -rf
  s:
    doc: |-
      -s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500
    type: string?
    inputBinding:
      prefix: -s=
      separate: false
  savd:
    doc: -savd=<rdatafile> OR -savd, save Rdata file
    type: boolean?
    inputBinding:
      valueFrom: |-
        ${ if (self) return "-savd=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.Rdata"; return null}
  savn:
    doc: -savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)
    type: boolean?
    inputBinding:
      valueFrom: |-
        ${ if (self) return "-savn=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.narrowPeak"; return null}
  savp:
    doc: save cross-correlation plot
    type: boolean?
    inputBinding:
      valueFrom: |-
        ${ if (self) return "-savp=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.pdf"; return null}
  savr:
    doc: |-
      -savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment)
    type: boolean?
    inputBinding:
      valueFrom: |-
        ${ if (self) return "-savr=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.regionPeak"; return null}
  speak:
    doc: -speak=<strPeak>, user-defined cross-correlation peak strandshift
    type: string?
    inputBinding:
      prefix: -speak=
      separate: false
  x:
    doc: |-
      -x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)
    type: string?
    inputBinding:
      prefix: -x=
      separate: false

outputs:
  output_spp_cross_corr:
    doc: peakshift/phantomPeak results summary file
    type: File
    outputBinding:
      glob: |-
        $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.txt")
  output_spp_cross_corr_plot:
    doc: peakshift/phantomPeak results summary plot
    type: File?
    outputBinding:
      glob: |-
        $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.pdf")
  output_spp_narrow_peak:
    doc: narrowPeak output file
    type: File?
    outputBinding:
      glob: |-
        $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.narrowPeak")
  output_spp_rdata:
    doc: Rdata file from the run_spp.R run
    type: File?
    outputBinding:
      glob: |-
        $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.Rdata")
  output_spp_region_peak:
    doc: regionPeak output file
    type: File?
    outputBinding:
      glob: |-
        $(inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp.regionPeak")

baseCommand: run_spp.R
arguments:
- valueFrom: |-
    $("-out=" + inputs.input_bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".spp_cross_corr.txt")
- valueFrom: $("-tmpdir="+runtime.tmpdir)
  shellQuote: false

hints:
  DockerRequirement:
    dockerPull: dukegcb/spp
