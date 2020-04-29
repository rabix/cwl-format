class: Workflow
cwlVersion: v1.0
doc: ChIP-seq - map - PCR Bottleneck Coefficients

requirements:
- class: ScatterFeatureRequirement

inputs:
  genome_sizes:
    type: File
  input_bam_files:
    type: File[]
  input_output_filenames:
    type: string[]

outputs:
  pbc_file:
    type: File[]
    outputSource: compute_pbc/pbc

steps:
  bedtools_genomecov:
    in:
      bg:
        default: true
      g: genome_sizes
      ibam: input_bam_files
    scatter: ibam
    run: execute_pcr_bottleneck_coef.cwl.steps/bedtools_genomecov.cwl
    out:
    - output_bedfile
  compute_pbc:
    in:
      bedgraph_file: bedtools_genomecov/output_bedfile
      output_filename: input_output_filenames
    scatter:
    - bedgraph_file
    - output_filename
    scatterMethod: dotproduct
    run: execute_pcr_bottleneck_coef.cwl.steps/compute_pbc.cwl
    out:
    - pbc
