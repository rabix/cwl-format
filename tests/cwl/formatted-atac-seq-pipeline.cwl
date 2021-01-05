#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: ATAC-seq-pipeline-se
doc: 'ATAC-seq pipeline - reads: SE'
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement

inputs:
  as_narrowPeak_file:
    doc: Definition narrowPeak file in AutoSql format (used in bedToBigBed)
    type: File
  default_adapters_file:
    doc: Adapters file
    type: File
  genome_effective_size:
    doc: |-
      Effective genome size used by MACS2. It can be numeric or a shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
    type: string
    default: hs
  genome_ref_first_index_file:
    doc: |-
      "First index file of Bowtie reference genome with extension 1.ebwt. \ (Note: the rest of the index files MUST be in the same folder)" 
    type: File
    secondaryFiles:
    - ^^.2.ebwt
    - ^^.3.ebwt
    - ^^.4.ebwt
    - ^^.rev.1.ebwt
    - ^^.rev.2.ebwt
  genome_sizes_file:
    doc: Genome sizes tab-delimited file (used in samtools)
    type: File
  input_fastq_files:
    type: File[]
  nthreads_map:
    doc: Number of threads required for the 03-map step
    type: int
  nthreads_peakcall:
    doc: Number of threads required for the 04-peakcall step
    type: int
  nthreads_qc:
    doc: Number of threads required for the 01-qc step
    type: int
  nthreads_quant:
    doc: Number of threads required for the 05-quantification step
    type: int
  nthreads_trimm:
    doc: Number of threads required for the 02-trim step
    type: int
  picard_jar_path:
    doc: Picard Java jar file
    type: string
  picard_java_opts:
    doc: |-
      JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
    type: string?
  trimmomatic_jar_path:
    doc: Trimmomatic Java jar file
    type: string
  trimmomatic_java_opts:
    doc: |-
      JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
    type: string?

outputs:
  map_bowtie_log_files:
    doc: Bowtie log file with mapping stats
    type: File[]
    outputSource: map/output_bowtie_log
  map_dedup_bam_files:
    doc: Filtered BAM files (post-processing end point)
    type: File[]
    outputSource: map/output_data_sorted_dups_marked_bam_files
  map_mark_duplicates_files:
    doc: |-
      Summary of duplicates removed with Picard tool MarkDuplicates (for multiple reads aligned to the same positions
    type: File[]
    outputSource: map/output_picard_mark_duplicates_files
  map_pbc_files:
    doc: PCR Bottleneck Coefficient files (used to flag samples when pbc<0.5)
    type: File[]
    outputSource: map/output_pbc_files
  map_percent_mitochondrial_reads:
    doc: Percentage of mitochondrial reads
    type: File[]
    outputSource: map/output_percent_mitochondrial_reads
  map_preseq_c_curve_files:
    doc: Preseq c_curve output files
    type: File[]
    outputSource: map/output_preseq_c_curve_files
  map_preseq_percentage_uniq_reads:
    doc: Preseq percentage of uniq reads
    type: File[]
    outputSource: map/output_percentage_uniq_reads
  map_read_count_mapped:
    doc: Read counts of the mapped BAM files
    type: File[]
    outputSource: map/output_read_count_mapped
  peakcall_extended_peak_file:
    doc: Extended fragment peaks in ENCODE Peak file format
    type: File[]
    outputSource: peak_call/output_extended_peak_file
  peakcall_filtered_read_count_file:
    doc: Filtered read count after peak calling
    type: File[]
    outputSource: peak_call/output_filtered_read_count_file
  peakcall_peak_bigbed_file:
    doc: Peaks in bigBed format
    type: File[]
    outputSource: peak_call/output_peak_bigbed_file
  peakcall_peak_count_within_replicate:
    doc: Peak counts within replicate
    type: File[]
    outputSource: peak_call/output_peak_count_within_replicate
  peakcall_peak_file:
    doc: Peaks in ENCODE Peak file format
    type: File[]
    outputSource: peak_call/output_peak_file
  peakcall_peak_summits_file:
    doc: Peaks summits in bedfile format
    type: File[]
    outputSource: peak_call/output_peak_summits_file
  peakcall_peak_xls_file:
    doc: Peak calling report file
    type: File[]
    outputSource: peak_call/output_peak_xls_file
  peakcall_read_in_peak_count_within_replicate:
    doc: Peak counts within replicate
    type: File[]
    outputSource: peak_call/output_read_in_peak_count_within_replicate
  peakcall_spp_x_cross_corr:
    doc: SPP strand cross correlation summary
    type: File[]
    outputSource: peak_call/output_spp_x_cross_corr
  peakcall_spp_x_cross_corr_plot:
    doc: SPP strand cross correlation plot
    type: File[]
    outputSource: peak_call/output_spp_cross_corr_plot
  qc_count_raw_reads:
    doc: Raw read counts of fastq files after QC
    type: File[]
    outputSource: qc/output_count_raw_reads
  qc_diff_counts:
    doc: Diff file between number of raw reads and number of reads counted by FASTQC,
    type: File[]
    outputSource: qc/output_diff_counts
  qc_fastqc_data_files:
    doc: FastQC data files
    type: File[]
    outputSource: qc/output_fastqc_data_files
  qc_fastqc_report_files:
    doc: FastQC reports in zip format
    type: File[]
    outputSource: qc/output_fastqc_report_files
  quant_bigwig_norm_files:
    doc: Normalized reads bigWig (signal) files
    type: File[]
    outputSource: quant/bigwig_norm_files
  quant_bigwig_raw_files:
    doc: Raw reads bigWig (signal) files
    type: File[]
    outputSource: quant/bigwig_raw_files
  trimm_fastq_files:
    doc: FASTQ files  after trimming
    type: File[]
    outputSource: trimm/output_data_fastq_trimmed_files
  trimm_raw_counts:
    doc: Raw read counts of fastq files after trimming
    type: File[]
    outputSource: trimm/output_trimmed_fastq_read_count

steps:
  map:
    in:
      genome_ref_first_index_file: genome_ref_first_index_file
      genome_sizes_file: genome_sizes_file
      input_fastq_files: trimm/output_data_fastq_trimmed_files
      nthreads: nthreads_map
      picard_jar_path: picard_jar_path
      picard_java_opts: picard_java_opts
    run:
      cwlVersion: v1.0
      class: Workflow
      doc: 'ATAC-seq 03 mapping - reads: SE'

      requirements:
      - class: ScatterFeatureRequirement
      - class: SubworkflowFeatureRequirement
      - class: StepInputExpressionRequirement
      - class: InlineJavascriptRequirement

      inputs:
        genome_ref_first_index_file:
          doc: |-
            Bowtie first index files for reference genome (e.g. *1.ebwt). The rest of the files should be in the same folder.
          type: File
          secondaryFiles:
          - ^^.2.ebwt
          - ^^.3.ebwt
          - ^^.4.ebwt
          - ^^.rev.1.ebwt
          - ^^.rev.2.ebwt
        genome_sizes_file:
          doc: Genome sizes tab-delimited file (used in samtools)
          type: File
        input_fastq_files:
          doc: Input fastq files
          type: File[]
        nthreads:
          type: int
          default: 1
        picard_jar_path:
          doc: Picard Java jar file
          type: string
          default: /usr/picard/picard.jar
        picard_java_opts:
          doc: |-
            JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
          type: string?

      outputs:
        output_bowtie_log:
          doc: Bowtie log file.
          type: File[]
          outputSource: bowtie-se/output_bowtie_log
        output_data_sorted_dedup_bam_files:
          doc: BAM files without duplicate reads.
          type: File[]
          outputSource: index_dedup_bams/indexed_file
        output_data_sorted_dups_marked_bam_files:
          doc: BAM files with marked duplicate reads.
          type: File[]
          outputSource: index_dups_marked_bams/indexed_file
        output_pbc_files:
          doc: PCR Bottleneck Coeficient files.
          type: File[]
          outputSource: execute_pcr_bottleneck_coef/pbc_file
        output_percent_mitochondrial_reads:
          doc: Percentage of mitochondrial reads.
          type: File[]
          outputSource: percent_mitochondrial_reads/percent_map
        output_percentage_uniq_reads:
          doc: Percentage of uniq reads from preseq c_curve output
          type: File[]
          outputSource: percent_uniq_reads/output
        output_picard_mark_duplicates_files:
          doc: Picard MarkDuplicates metrics files.
          type: File[]
          outputSource: mark_duplicates/output_metrics_file
        output_preseq_c_curve_files:
          doc: Preseq c_curve output files.
          type: File[]
          outputSource: preseq-c-curve/output_file
        output_read_count_mapped:
          doc: Read counts of the mapped BAM files
          type: File[]
          outputSource: mapped_reads_count/output
        output_read_count_mapped_filtered:
          doc: Read counts of the mapped and filtered BAM files
          type: File[]
          outputSource: mapped_filtered_reads_count/output_read_count

      steps:
        bam_idxstats:
          in:
            bam: index_bams/indexed_file
          scatter: bam
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
          out:
          - idxstats_file
        bowtie-se:
          in:
            X:
              valueFrom: ${return 2000}
            genome_ref_first_index_file: genome_ref_first_index_file
            input_fastq_file: input_fastq_files
            nthreads: nthreads
            output_filename: extract_basename_2/output_path
            v:
              valueFrom: ${return 2}
          scatter:
          - input_fastq_file
          - output_filename
          scatterMethod: dotproduct
          run:
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
          out:
          - output_aligned_file
          - output_bowtie_log
        execute_pcr_bottleneck_coef:
          in:
            genome_sizes: genome_sizes_file
            input_bam_files: filtered2sorted/sorted_file
            input_output_filenames: extract_basename_2/output_path
          run:
            cwlVersion: v1.0
            class: Workflow
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
                run:
                  cwlVersion: v1.0
                  class: CommandLineTool
                  doc: |-
                    Tool:    bedtools genomecov (aka genomeCoverageBed)
                    Version: v2.25.0
                    Summary: Compute the coverage of a feature file among a genome.

                    Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

                    Options: 
                    	-ibam		The input file is in BAM format.
                    			Note: BAM _must_ be sorted by position

                    	-d		Report the depth at each genome position (with one-based coordinates).
                    			Default behavior is to report a histogram.

                    	-dz		Report the depth at each genome position (with zero-based coordinates).
                    			Reports only non-zero positions.
                    			Default behavior is to report a histogram.

                    	-bg		Report depth in BedGraph format. For details, see:
                    			genome.ucsc.edu/goldenPath/help/bedgraph.html

                    	-bga		Report depth in BedGraph format, as above (-bg).
                    			However with this option, regions with zero 
                    			coverage are also reported. This allows one to
                    			quickly extract all regions of a genome with 0 
                    			coverage by applying: "grep -w 0$" to the output.

                    	-split		Treat "split" BAM or BED12 entries as distinct BED intervals.
                    			when computing coverage.
                    			For BAM files, this uses the CIGAR "N" and "D" operations 
                    			to infer the blocks for computing coverage.
                    			For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
                    			fields (i.e., columns 10,11,12).

                    	-strand		Calculate coverage of intervals from a specific strand.
                    			With BED files, requires at least 6 columns (strand is column 6). 
                    			- (STRING): can be + or -

                    	-5		Calculate coverage of 5" positions (instead of entire interval).

                    	-3		Calculate coverage of 3" positions (instead of entire interval).

                    	-max		Combine all positions with a depth >= max into
                    			a single bin in the histogram. Irrelevant
                    			for -d and -bedGraph
                    			- (INTEGER)

                    	-scale		Scale the coverage by a constant factor.
                    			Each coverage value is multiplied by this factor before being reported.
                    			Useful for normalizing coverage by, e.g., reads per million (RPM).
                    			- Default is 1.0; i.e., unscaled.
                    			- (FLOAT)

                    	-trackline	Adds a UCSC/Genome-Browser track line definition in the first line of the output.
                    			- See here for more details about track line definition:
                    			      http://genome.ucsc.edu/goldenPath/help/bedgraph.html
                    			- NOTE: When adding a trackline definition, the output BedGraph can be easily
                    			      uploaded to the Genome Browser as a custom track,
                    			      BUT CAN NOT be converted into a BigWig file (w/o removing the first line).

                    	-trackopts	Writes additional track line definition parameters in the first line.
                    			- Example:
                    			   -trackopts 'name="My Track" visibility=2 color=255,30,30'
                    			   Note the use of single-quotes if you have spaces in your parameters.
                    			- (TEXT)

                    Notes: 
                    	(1) The genome file should tab delimited and structured as follows:
                    	 <chromName><TAB><chromSize>

                    	For example, Human (hg19):
                    	chr1	249250621
                    	chr2	243199373
                    	...
                    	chr18_gl000207_random	4262

                    	(2) The input BED (-i) file must be grouped by chromosome.
                    	 A simple "sort -k 1,1 <BED> > <BED>.sorted" will suffice.

                    	(3) The input BAM (-ibam) file must be sorted by position.
                    	 A "samtools sort <BAM>" should suffice.

                    Tips: 
                    	One can use the UCSC Genome Browser's MySQL database to extract
                    	chromosome sizes. For example, H. sapiens:

                    	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                    	"select chrom, size from hg19.chromInfo" > hg19.genome

                  requirements:
                    InlineJavascriptRequirement: {}
                    ShellCommandRequirement: {}

                  inputs:
                    '3':
                      doc: "\tCalculate coverage of 3\" positions (instead of entire\
                        \ interval).\n"
                      type: boolean?
                      inputBinding:
                        prefix: '-3'
                        position: 1
                    '5':
                      doc: "\tCalculate coverage of 5\" positions (instead of entire\
                        \ interval).\n"
                      type: boolean?
                      inputBinding:
                        prefix: '-5'
                        position: 1
                    bg:
                      doc: |
                        	Report depth in BedGraph format. For details, see:
                        genome.ucsc.edu/goldenPath/help/bedgraph.html
                      type: boolean?
                      inputBinding:
                        prefix: -bg
                        position: 1
                    bga:
                      doc: |
                        	Report depth in BedGraph format, as above (-bg).
                        However with this option, regions with zero
                        coverage are also reported. This allows one to
                        quickly extract all regions of a genome with 0
                        coverage by applying: "grep -w 0$" to the output.
                      type: boolean?
                      inputBinding:
                        prefix: -bga
                        position: 1
                    d:
                      doc: |
                        	Report the depth at each genome position (with one-based coordinates).
                        Default behavior is to report a histogram.
                      type: boolean?
                      inputBinding:
                        prefix: -d
                        position: 1
                    dz:
                      doc: |
                        	Report the depth at each genome position (with zero-based coordinates).
                        Reports only non-zero positions.
                        Default behavior is to report a histogram.
                      type: boolean?
                      inputBinding:
                        prefix: -dz
                        position: 1
                    g:
                      doc: <genome sizes>
                      type: File
                      inputBinding:
                        prefix: -g
                        position: 3
                    ibam:
                      doc: "\tThe input file is in BAM format.\nNote: BAM _must_ be\
                        \ sorted by position\n"
                      type: File
                      inputBinding:
                        prefix: -ibam
                        position: 2
                    max:
                      doc: |
                        	Combine all positions with a depth >= max into
                        a single bin in the histogram. Irrelevant
                        for -d and -bedGraph
                        - (INTEGER)
                      type: int?
                      inputBinding:
                        prefix: -max
                        position: 1
                    scale:
                      doc: |
                        	Scale the coverage by a constant factor.
                        Each coverage value is multiplied by this factor before being reported.
                        Useful for normalizing coverage by, e.g., reads per million (RPM).
                        - Default is 1.0; i.e., unscaled.
                        - (FLOAT)
                      type: float?
                      inputBinding:
                        prefix: -scale
                        position: 1
                    split:
                      doc: |
                        	Treat "split" BAM or BED12 entries as distinct BED intervals.
                        when computing coverage.
                        For BAM files, this uses the CIGAR "N" and "D" operations
                        to infer the blocks for computing coverage.
                        For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
                        fields (i.e., columns 10,11,12).
                      type: boolean?
                      inputBinding:
                        prefix: -split
                        position: 1
                    strand:
                      doc: |
                        	Calculate coverage of intervals from a specific strand.
                        With BED files, requires at least 6 columns (strand is column 6).
                        - (STRING): can be + or -
                      type: string?
                      inputBinding:
                        prefix: -strand
                        position: 1
                    trackline:
                      doc: |
                        Adds a UCSC/Genome-Browser track line definition in the first line of the output.
                        - See here for more details about track line definition:
                        http://genome.ucsc.edu/goldenPath/help/bedgraph.html
                        - NOTE: When adding a trackline definition, the output BedGraph can be easily
                        uploaded to the Genome Browser as a custom track,
                        BUT CAN NOT be converted into a BigWig file (w/o removing the first line).
                      type: boolean?
                      inputBinding:
                        prefix: -trackline
                        position: 1
                    trackopts:
                      doc: |
                        Writes additional track line definition parameters in the first line.
                        - Example:
                        -trackopts 'name="My Track" visibility=2 color=255,30,30'
                        Note the use of single-quotes if you have spaces in your parameters.
                        - (TEXT)
                      type: string?
                      inputBinding:
                        prefix: -trackopts
                        position: 1

                  outputs:
                    output_bedfile:
                      type: File
                      outputBinding:
                        glob: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
                          '') + '.bdg')
                  stdout: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
                    '') + '.bdg')

                  baseCommand:
                  - bedtools
                  - genomecov

                  hints:
                    DockerRequirement:
                      dockerPull: dukegcb/bedtools
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
                run:
                  cwlVersion: v1.0
                  class: CommandLineTool
                  doc: Compute PCR Bottleneck Coeficient from BedGraph file.

                  inputs:
                    bedgraph_file:
                      type: File
                      inputBinding:
                        position: 1
                    output_filename:
                      type: string

                  outputs:
                    pbc:
                      type: File
                      outputBinding:
                        glob: $(inputs.output_filename + '.PBC.txt')
                  stdout: $(inputs.output_filename + '.PBC.txt')

                  baseCommand:
                  - awk
                  - $4==1 {N1 += $3 - $2}; $4>=1 {Nd += $3 - $2} END {print N1/Nd}
                out:
                - pbc
          out:
          - pbc_file
        extract_basename_1:
          in:
            input_file: input_fastq_files
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts the base name of a file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_basename:
                type: string
                outputBinding:
                  outputEval: |-
                    $(inputs.input_file.path.substr(inputs.input_file.path.lastIndexOf('/') + 1, inputs.input_file.path.lastIndexOf('.') - (inputs.input_file.path.lastIndexOf('/') + 1)))

            baseCommand: echo

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_basename
        extract_basename_2:
          in:
            file_path: extract_basename_1/output_basename
          scatter: file_path
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts the base name of a file

            inputs:
              file_path:
                type: string
                inputBinding:
                  position: 1

            outputs:
              output_path:
                type: string
                outputBinding:
                  outputEval: $(inputs.file_path.replace(/\.[^/.]+$/, ""))

            baseCommand: echo

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_path
        filter-unmapped:
          in:
            input_file: sort_bams/sorted_file
            output_filename: extract_basename_2/output_path
          scatter:
          - input_file
          - output_filename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                doc: Aligned file to be sorted with samtools
                type: File
                inputBinding:
                  position: 1000
              nthreads:
                doc: Number of threads used in sorting
                type: int
                default: 1
                inputBinding:
                  prefix: -@
                  position: 1
              output_filename:
                doc: Basename for the output file
                type: string

            outputs:
              filtered_file:
                doc: Filter unmapped reads in aligned file
                type: File
                outputBinding:
                  glob: $(inputs.output_filename + '.accepted_hits.bam')
            stdout: $(inputs.output_filename + '.accepted_hits.bam')

            baseCommand:
            - samtools
            - view
            - -F
            - '4'
            - -b
            - -h

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - filtered_file
        filtered2sorted:
          in:
            input_file: filter-unmapped/filtered_file
            nthreads: nthreads
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: suffix of the transformed SAM/BAM file (including extension,
                  e.g. .filtered.sam)
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
          out:
          - sorted_file
        index_bams:
          in:
            input_file: sort_bams/sorted_file
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InitialWorkDirRequirement:
                listing:
                - $(inputs.input_file)
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                doc: Aligned file to be sorted with samtools
                type: File
                inputBinding:
                  position: 1

            outputs:
              indexed_file:
                doc: Indexed BAM file
                type: File
                secondaryFiles: .bai
                outputBinding:
                  glob: $(inputs.input_file.basename)

            baseCommand:
            - samtools
            - index
            arguments:
            - position: 2
              valueFrom: $(inputs.input_file.basename + '.bai')

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - indexed_file
        index_dedup_bams:
          in:
            input_file: sort_dedup_bams/sorted_file
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InitialWorkDirRequirement:
                listing:
                - $(inputs.input_file)
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                doc: Aligned file to be sorted with samtools
                type: File
                inputBinding:
                  position: 1

            outputs:
              indexed_file:
                doc: Indexed BAM file
                type: File
                secondaryFiles: .bai
                outputBinding:
                  glob: $(inputs.input_file.basename)

            baseCommand:
            - samtools
            - index
            arguments:
            - position: 2
              valueFrom: $(inputs.input_file.basename + '.bai')

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - indexed_file
        index_dups_marked_bams:
          in:
            input_file: sort_dups_marked_bams/sorted_file
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InitialWorkDirRequirement:
                listing:
                - $(inputs.input_file)
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                doc: Aligned file to be sorted with samtools
                type: File
                inputBinding:
                  position: 1

            outputs:
              indexed_file:
                doc: Indexed BAM file
                type: File
                secondaryFiles: .bai
                outputBinding:
                  glob: $(inputs.input_file.basename)

            baseCommand:
            - samtools
            - index
            arguments:
            - position: 2
              valueFrom: $(inputs.input_file.basename + '.bai')

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - indexed_file
        index_filtered_bam:
          in:
            input_file: filtered2sorted/sorted_file
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InitialWorkDirRequirement:
                listing:
                - $(inputs.input_file)
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                doc: Aligned file to be sorted with samtools
                type: File
                inputBinding:
                  position: 1

            outputs:
              indexed_file:
                doc: Indexed BAM file
                type: File
                secondaryFiles: .bai
                outputBinding:
                  glob: $(inputs.input_file.basename)

            baseCommand:
            - samtools
            - index
            arguments:
            - position: 2
              valueFrom: $(inputs.input_file.basename + '.bai')

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - indexed_file
        mapped_filtered_reads_count:
          in:
            input_bam_file: sort_dedup_bams/sorted_file
            output_suffix:
              valueFrom: .mapped_and_filtered.read_count.txt
          scatter: input_bam_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
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
          out:
          - output_read_count
        mapped_reads_count:
          in:
            bowtie_log: bowtie-se/output_bowtie_log
          scatter: bowtie_log
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Get number of processed reads from Bowtie log.

            requirements:
              InlineJavascriptRequirement: {}
              ShellCommandRequirement: {}

            inputs:
              bowtie_log:
                type: File
                inputBinding: {}

            outputs:
              output:
                type: File
                outputBinding:
                  glob: $(inputs.bowtie_log.path.replace(/^.*[\\\/]/, '') + '.read_count.mapped')
            stdout: $(inputs.bowtie_log.path.replace(/^.*[\\\/]/, '') + '.read_count.mapped')

            baseCommand: read-count-from-bowtie-log.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output
        mark_duplicates:
          in:
            input_file: index_filtered_bam/indexed_file
            java_opts: picard_java_opts
            output_filename: extract_basename_2/output_path
            output_suffix:
              valueFrom: bam
            picard_jar_path: picard_jar_path
          scatter:
          - input_file
          - output_filename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: One or more input SAM or BAM files to analyze. Must be coordinate
                  sorted.
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
          out:
          - output_metrics_file
          - output_dedup_bam_file
        percent_mitochondrial_reads:
          in:
            chrom:
              valueFrom: chrM
            idxstats: bam_idxstats/idxstats_file
            output_filename:
              valueFrom: |-
                ${return inputs.idxstats.basename.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '').replace(/\.[^/.]+$/, '').replace(/\.[^/.]+$/, '.mitochondrial_percentage.txt')}
          scatter: idxstats
          run:
            cwlVersion: v1.0
            class: ExpressionTool

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              chrom:
                doc: Query chromosome used to calculate percentage
                type: string
              idxstats:
                doc: Samtools idxstats file
                type: File
                inputBinding:
                  loadContents: true
              output_filename:
                doc: Save the percentage in a file of the given name
                type: string?

            outputs:
              percent_map:
                type:
                - File
                - string
            expression: |
              ${
                var regExp = new RegExp(inputs.chrom + "\\s\\d+\\s(\\d+)\\s(\\d+)");
                var match = inputs.idxstats.contents.match(regExp);
                if (match){
                  var chrom_mapped_reads = match[1];
                  var total_reads = inputs.idxstats.contents.split("\n")
                    .map(function(x){
                      var rr = x.match(/.*\s\d+\s(\d+)\s\d+/);
                      return (rr ? rr[1] : 0);
                    })
                    .reduce(function(a, b) { return Number(a) + Number(b); });

                  var output = (100*chrom_mapped_reads/total_reads).toFixed(4) + "%" + "\n";

                  if (inputs.output_filename){
                    return {
                      percent_map : {
                        "class": "File",
                        "basename" : inputs.output_filename,
                        "contents" : output,
                      }
                    }
                  }
                  return output;
                }
              }
          out:
          - percent_map
        percent_uniq_reads:
          in:
            preseq_c_curve_outfile: preseq-c-curve/output_file
          scatter: preseq_c_curve_outfile
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Get number of processed reads from Bowtie log.

            requirements:
              InlineJavascriptRequirement: {}
              ShellCommandRequirement: {}

            inputs:
              preseq_c_curve_outfile:
                type: File
                inputBinding: {}

            outputs:
              output:
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.preseq_c_curve_outfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + '.percentage_unique_reads.txt')
            stdout: |-
              $(inputs.preseq_c_curve_outfile.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + '.percentage_unique_reads.txt')

            baseCommand: percent-uniq-reads-from-preseq.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output
        preseq-c-curve:
          in:
            input_sorted_file: filtered2sorted/sorted_file
            output_file_basename: extract_basename_2/output_path
          scatter:
          - input_sorted_file
          - output_file_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |-
              Usage: c_curve [OPTIONS] <sorted-bed-file>

              Options:
                -o, -output   yield output file (default: stdout) 
                -s, -step     step size in extrapolations (default: 1e+06) 
                -v, -verbose  print more information 
                -P, -pe       input is paired end read file 
                -H, -hist     input is a text file containing the observed histogram 
                -V, -vals     input is a text file containing only the observed counts 
                -B, -bam      input is in BAM format 
                -l, -seg_len  maximum segment length when merging paired end bam reads 
                              (default: 5000) 

              Help options:
                -?, -help     print this help message 
                    -about    print about message

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              B:
                doc: "-bam      input is in BAM format \n"
                type: boolean
                default: true
                inputBinding:
                  prefix: -B
                  position: 1
              H:
                doc: "-hist     input is a text file containing the observed histogram\
                  \ \n"
                type: File?
                inputBinding:
                  prefix: -H
                  position: 1
              V:
                doc: "-vals     input is a text file containing only the observed\
                  \ counts \n"
                type: File?
                inputBinding:
                  prefix: -V
                  position: 1
              input_sorted_file:
                doc: Sorted bed or BAM file
                type: File
                inputBinding:
                  position: 2
              l:
                doc: |
                  -seg_len  maximum segment length when merging paired end bam reads 
                  (default: 5000)
                  Help options:
                  -?, -help     print this help message
                  -about    print about message
                type: int?
                inputBinding:
                  prefix: -l
                  position: 1
              output_file_basename:
                type: string
              pe:
                doc: "-pe       input is paired end read file \n"
                type: boolean?
                inputBinding:
                  prefix: -P
                  position: 1
              s:
                doc: "-step     step size in extrapolations (default: 1e+06) \n"
                type: float?
                inputBinding:
                  prefix: -s
                  position: 1
              v:
                doc: "-verbose  print more information \n"
                type: boolean
                default: false
                inputBinding:
                  prefix: -v
                  position: 1

            outputs:
              output_file:
                type: File
                outputBinding:
                  glob: $(inputs.output_file_basename + '.preseq_c_curve.txt')
            stdout: $(inputs.output_file_basename + '.preseq_c_curve.txt')

            baseCommand:
            - preseq
            - c_curve

            hints:
              DockerRequirement:
                dockerPull: reddylab/preseq:2.0
          out:
          - output_file
        remove_duplicates:
          in:
            F:
              valueFrom: ${return 1024}
            b:
              valueFrom: ${return true}
            input_file: index_dups_marked_bams/indexed_file
            outfile_name:
              valueFrom: ${return inputs.input_file.basename.replace('dups_marked',
                'dedup')}
            suffix:
              valueFrom: .dedup.bam
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              F:
                doc: only include reads with none of the bits set in INT set in FLAG
                  [0]
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
                doc: suffix of the transformed SAM/BAM file (including extension,
                  e.g. .filtered.sam)
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
          out:
          - outfile
        sam2bam:
          in:
            input_file: bowtie-se/output_aligned_file
            nthreads: nthreads
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              S:
                doc: Input format autodetected
                type: boolean
                default: true
                inputBinding:
                  prefix: -S
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

            outputs:
              bam_file:
                doc: Aligned file in BAM format
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.bam')
            stdout: |-
              $(inputs.input_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.bam')

            baseCommand:
            - samtools
            - view
            - -b

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - bam_file
        sort_bams:
          in:
            input_file: sam2bam/bam_file
            nthreads: nthreads
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: suffix of the transformed SAM/BAM file (including extension,
                  e.g. .filtered.sam)
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
          out:
          - sorted_file
        sort_dedup_bams:
          in:
            input_file: remove_duplicates/outfile
            nthreads: nthreads
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: suffix of the transformed SAM/BAM file (including extension,
                  e.g. .filtered.sam)
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
          out:
          - sorted_file
        sort_dups_marked_bams:
          in:
            input_file: mark_duplicates/output_dedup_bam_file
            nthreads: nthreads
            suffix:
              valueFrom: .dups_marked.bam
          scatter:
          - input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: suffix of the transformed SAM/BAM file (including extension,
                  e.g. .filtered.sam)
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
          out:
          - sorted_file
    out:
    - output_data_sorted_dedup_bam_files
    - output_data_sorted_dups_marked_bam_files
    - output_picard_mark_duplicates_files
    - output_pbc_files
    - output_bowtie_log
    - output_preseq_c_curve_files
    - output_percentage_uniq_reads
    - output_read_count_mapped
    - output_percent_mitochondrial_reads
  peak_call:
    in:
      as_narrowPeak_file: as_narrowPeak_file
      genome_effective_size: genome_effective_size
      input_bam_files: map/output_data_sorted_dedup_bam_files
      input_bam_format:
        valueFrom: BAM
      input_genome_sizes: genome_sizes_file
      nthreads: nthreads_peakcall
    run:
      cwlVersion: v1.0
      class: Workflow
      doc: ATAC-seq 04 quantification - SE

      requirements:
      - class: ScatterFeatureRequirement
      - class: StepInputExpressionRequirement
      - class: InlineJavascriptRequirement

      inputs:
        as_narrowPeak_file:
          doc: Definition narrowPeak file in AutoSql format (used in bedToBigBed)
          type: File
        genome_effective_size:
          doc: |-
            Effective genome size used by MACS2. It can be numeric or a shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8), Default:hs
          type: string
          default: hs
        input_bam_files:
          type: File[]
        input_genome_sizes:
          doc: Two column tab-delimited file with chromosome size information
          type: File
        nthreads:
          type: int
          default: 1

      outputs:
        output_extended_peak_file:
          doc: peakshift/phantomPeak extended fragment results file
          type: File[]
          outputSource: peak-calling/output_ext_frag_bdg_file
        output_filtered_read_count_file:
          doc: Filtered read count reported by MACS2
          type: File[]
          outputSource: count-reads-filtered/read_count_file
        output_peak_bigbed_file:
          doc: Peaks in bigBed format
          type: File[]
          outputSource: peaks-bed-to-bigbed/bigbed
        output_peak_count_within_replicate:
          doc: Peak counts within replicate
          type: File[]
          outputSource: count-peaks/output_counts
        output_peak_file:
          doc: peakshift/phantomPeak results file
          type: File[]
          outputSource: peak-calling/output_peak_file
        output_peak_summits_file:
          doc: File containing peak summits
          type: File[]
          outputSource: peak-calling/output_peak_summits_file
        output_peak_xls_file:
          doc: Peak calling report file (*_peaks.xls file produced by MACS2)
          type: File[]
          outputSource: peak-calling/output_peak_xls_file
        output_read_in_peak_count_within_replicate:
          doc: Reads peak counts within replicate
          type: File[]
          outputSource: extract-count-reads-in-peaks/output_read_count
        output_spp_cross_corr_plot:
          doc: peakshift/phantomPeak results file
          type: File[]
          outputSource: spp/output_spp_cross_corr_plot
        output_spp_x_cross_corr:
          doc: peakshift/phantomPeak results file
          type: File[]
          outputSource: spp/output_spp_cross_corr

      steps:
        count-peaks:
          in:
            input_file: peak-calling/output_peak_file
            output_suffix:
              valueFrom: .peak_count.within_replicate.txt
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Counts lines in a file and returns a suffixed file with that number

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                type: File
              output_suffix:
                type: string
                default: .count

            outputs:
              output_counts:
                type: File
                outputBinding:
                  glob: $(inputs.input_file.path.replace(/^.*[\\\/]/, '') + inputs.output_suffix)
            stdout: $(inputs.input_file.path.replace(/^.*[\\\/]/, '') + inputs.output_suffix)

            baseCommand:
            - wc
            - -l
            stdin: $(inputs.input_file.path)
          out:
          - output_counts
        count-reads-filtered:
          in:
            peak_xls_file: peak-calling/output_peak_xls_file
          scatter: peak_xls_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Count number of dedup-ed reads used in peak calling

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              peak_xls_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              read_count_file:
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.peak_xls_file.path.replace(/^.*[\\\/]/, '').replace(/\_peaks\.xls$/, '_read_count.txt'))
            stdout: |-
              $(inputs.peak_xls_file.path.replace(/^.*[\\\/]/, '').replace(/\_peaks\.xls$/, '_read_count.txt'))

            baseCommand: count-filtered-reads-macs2.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - read_count_file
        extract-count-reads-in-peaks:
          in:
            input_bam_file: filter-reads-in-peaks/filtered_file
            output_suffix:
              valueFrom: .read_count.within_replicate.txt
          scatter: input_bam_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
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
          out:
          - output_read_count
        extract-peak-frag-length:
          in:
            input_spp_txt_file: spp/output_spp_cross_corr
          scatter: input_spp_txt_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts best fragment length from SPP output text file

            inputs:
              input_spp_txt_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_best_frag_length:
                type: float
                outputBinding:
                  glob: best_frag_length
                  outputEval: $(Number(self[0].contents.replace('\n', '')))
                  loadContents: true
            stdout: best_frag_length

            baseCommand: extract-best-frag-length.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_best_frag_length
        filter-reads-in-peaks:
          in:
            input_bam_file: input_bam_files
            input_bedfile: peak-calling/output_peak_file
          scatter:
          - input_bam_file
          - input_bedfile
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Filter BAM file to only include reads overlapping with a BED file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_bam_file:
                doc: Aligned BAM file to filter
                type: File
                inputBinding:
                  position: 3
              input_bedfile:
                doc: Bedfile used to only include reads overlapping this BED FILE
                type: File
                inputBinding:
                  prefix: -L
                  position: 2

            outputs:
              filtered_file:
                doc: Filtered aligned BAM file by BED coordinates file
                type: File
                outputBinding:
                  glob: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '') + '.in_peaks.bam')
            stdout: $(inputs.input_bam_file.path.replace(/^.*[\\\/]/, '') + '.in_peaks.bam')

            baseCommand:
            - samtools
            - view
            - -b
            - -h

            hints:
              DockerRequirement:
                dockerPull: dukegcb/samtools:1.3
          out:
          - filtered_file
        peak-calling:
          in:
            format:
              valueFrom: BAM
            bdg:
              valueFrom: ${return true}
            extsize:
              valueFrom: ${return 200}
            g: genome_effective_size
            nomodel:
              valueFrom: ${return true}
            q:
              valueFrom: ${return 0.1}
            shift:
              valueFrom: ${return -100}
            treatment:
              valueFrom: $([self])
              source: input_bam_files
          scatter:
          - treatment
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              format:
                doc: |-
                  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file, "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or "BAMPE". The default AUTO option will let MACS decide which format the file is. Note that MACS can't detect "BAMPE" or "BEDPE" format with "AUTO", and you have to implicitly specify the format for "BAMPE" and "BEDPE". DEFAULT: "AUTO".
                type: string?
                inputBinding:
                  prefix: -f
                  position: 1
              SPMR:
                doc: |-
                  If True, MACS will save signal per million reads for fragment pileup profiles. Require --bdg to be set. Default: False 
                type: boolean?
                inputBinding:
                  prefix: --SPMR
                  position: 1
              bdg:
                doc: |-
                  Whether or not to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file. DEFAULT: True
                type: boolean?
                inputBinding:
                  prefix: --bdg
                  position: 1
              broad:
                doc: |-
                  If set, MACS will try to call broad peaks by linking nearby highly enriched regions. The linking region is controlled by another cutoff through --linking-cutoff. The maximum linking region length is 4 times of d from MACS. DEFAULT: False 
                type: boolean?
                inputBinding:
                  prefix: --broad
                  position: 1
              broad-cutoff:
                doc: |
                  BROADCUTOFF
                  Cutoff for broad region. This option is not available
                  unless --broad is set. If -p is set, this is a pvalue
                  cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1
                type: float?
                inputBinding:
                  prefix: --broad-cutoff
                  position: 1
              buffer-size:
                doc: |
                  BUFFER_SIZE
                  Buffer size for incrementally increasing internal
                  array size to store reads alignment information. In
                  most cases, you don't have to change this parameter.
                  However, if there are large number of
                  chromosomes/contigs/scaffolds in your alignment, it's
                  recommended to specify a smaller buffer size in order
                  to decrease memory usage (but it will take longer time
                  to read alignment files). Minimum memory requested for
                  reading an alignment file is about # of CHROMOSOME *
                  BUFFER_SIZE * 2 Bytes. DEFAULT: 100000
                type: int?
                inputBinding:
                  prefix: --buffer-size
                  position: 1
              bw:
                doc: |
                  BW               Band width for picking regions to compute fragment
                  size. This value is only used while building the
                  shifting model. DEFAULT: 300
                type: int?
                inputBinding:
                  prefix: --bw
                  position: 1
              call-summits:
                doc: |-
                  If set, MACS will use a more sophisticated signal processing approach to find subpeak summits in each enriched peak region. DEFAULT: False 
                type: boolean?
                inputBinding:
                  prefix: --call-summits
                  position: 1
              control:
                doc: Control sample file.
                type: File?
                inputBinding:
                  prefix: --control
                  position: 2
              cutoff-analysis:
                doc: |-
                  While set, MACS2 will analyze number or total length of peaks that can be called by different p-value cutoff then output a summary table to help user decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file. Note, minlen and maxgap may affect the results. WARNING: May take ~30 folds longer time to finish. DEFAULT: False Post-processing options: 
                type: boolean?
                inputBinding:
                  prefix: --cutoff-analysis
                  position: 1
              down-sample:
                doc: |-
                  When set, random sampling method will scale down the bigger sample. By default, MACS uses linear scaling. Warning: This option will make your result unstable and irreproducible since each time, random reads would be selected. Consider to use 'randsample' script instead. <not implmented>If used together with --SPMR, 1 million unique reads will be randomly picked.</not implemented> Caution: due to the implementation, the final number of selected reads may not be as you expected! DEFAULT: False 
                type: boolean?
                inputBinding:
                  prefix: --down-sample
                  position: 1
              extsize:
                doc: |-
                  The arbitrary extension size in bp. When nomodel is  true, MACS will use this value as fragment size to  extend each read towards 3' end, then pile them up.  It's exactly twice the number of obsolete SHIFTSIZE.  In previous language, each read is moved 5'->3'  direction to middle of fragment by 1/2 d, then  extended to both direction with 1/2 d. This is  equivalent to say each read is extended towards 5'->3'  into a d size fragment. DEFAULT: 200. EXTSIZE and  SHIFT can be combined when necessary. Check SHIFT  option.
                type: float?
                inputBinding:
                  prefix: --extsize
                  position: 1
              fe-cutoff:
                doc: |
                  FECUTOFF  When set, the value will be used to filter out peaks
                  with low fold-enrichment. Note, MACS2 use 1.0 as
                  pseudocount while calculating fold-enrichment.
                  DEFAULT: 1.0
                type: float?
                inputBinding:
                  prefix: --fe-cutoff
                  position: 1
              fix-bimodal:
                doc: |-
                  Whether turn on the auto pair model process. If set, when MACS failed to build paired model, it will use the nomodel settings, the --exsize parameter to extend each tags towards 3' direction. Not to use this automate fixation is a default behavior now. DEFAULT: False 
                type: boolean?
                inputBinding:
                  prefix: --fix-bimodal
                  position: 1
              g:
                doc: |-
                  Effective genome size. It can be 1.0e+9 or 1000000000,  or shortcuts:'hs' for human (2.7e9), 'mm' for mouse  (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for  fruitfly (1.2e8), Default:hs.
                type: string?
                inputBinding:
                  prefix: -g
                  position: 1
              keep-dup:
                doc: |
                  KEEPDUPLICATES
                  It controls the MACS behavior towards duplicate tags
                  at the exact same location -- the same coordination
                  and the same strand. The 'auto' option makes MACS
                  calculate the maximum tags at the exact same location
                  based on binomal distribution using 1e-5 as pvalue
                  cutoff; and the 'all' option keeps every tags. If an
                  integer is given, at most this number of tags will be
                  kept at the same location. The default is to keep one
                  tag at the same location. Default: 1
                type: string?
                inputBinding:
                  prefix: --keep-dup
                  position: 1
              llocal:
                doc: |
                  LARGELOCAL   The large nearby region in basepairs to calculate
                  dynamic lambda. This is used to capture the surround
                  bias. If you set this to 0, MACS will skip llocal
                  lambda calculation. *Note* that MACS will always
                  perform a d-size local lambda calculation. The final
                  local bias should be the maximum of the lambda value
                  from d, slocal, and llocal size windows. DEFAULT:
                  10000.
                type: int?
                inputBinding:
                  prefix: --llocal
                  position: 1
              m:
                doc: |
                  MFOLD MFOLD, --mfold MFOLD MFOLD
                  Select the regions within MFOLD range of high-
                  confidence enrichment ratio against background to
                  build model. Fold-enrichment in regions must be lower
                  than upper limit, and higher than the lower limit. Use
                  as "-m 10 30". DEFAULT:5 50
                type: string?
                inputBinding:
                  prefix: -m
                  position: 1
              nolambda:
                doc: |-
                  If True, MACS will use fixed background lambda as local lambda for every peak region. Normally, MACS calculates a dynamic local lambda to reflect the local bias due to potential chromatin structure. 
                type: boolean?
                inputBinding:
                  prefix: --nolambda
                  position: 1
              nomodel:
                doc: |-
                  	 Whether or not to build the shifting model. If True,  MACS will not build model. by default it means  shifting size = 100, try to set extsize to change it.  DEFAULT: False
                type: boolean?
                inputBinding:
                  prefix: --nomodel
                  position: 1
              p:
                doc: |-
                  Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually exclusive. If pvalue cutoff is  set, qvalue will not be calculated and reported as -1  in the final .xls file..
                type: float?
                inputBinding:
                  prefix: -p
                  position: 1
              q:
                doc: |-
                  Minimum FDR (q-value) cutoff for peak detection. DEFAULT: 0.05. -q, and -p are mutually exclusive.
                type: float?
                inputBinding:
                  prefix: -q
                  position: 1
              ratio:
                doc: |
                  RATIO         When set, use a custom scaling ratio of ChIP/control
                  (e.g. calculated using NCIS) for linear scaling.
                  DEFAULT: ingore
                type: float?
                inputBinding:
                  prefix: --ratio
                  position: 1
              s:
                doc: |
                  TSIZE, --tsize TSIZE
                  Tag size. This will overide the auto detected tag
                  size. DEFAULT: Not set
                type: int?
                inputBinding:
                  prefix: -s
                  position: 1
              seed:
                doc: |
                  SEED           Set the random seed while down sampling data. Must be
                  a non-negative integer in order to be effective.
                  DEFAULT: not set
                type: int?
                inputBinding:
                  prefix: --seed
                  position: 1
              shift:
                doc: |-
                  (NOT the legacy --shiftsize option!) The arbitrary shift in bp. Use discretion while setting it other than default value. When NOMODEL is set, MACS will use this value to move cutting ends (5') towards 5'->3' direction then apply EXTSIZE to extend them to fragments. When this value is negative, ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for ChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for detecting enriched cutting loci such as certain DNAseI-Seq datasets. Note, you can't set values other than 0 if format is BAMPE for paired-end data. DEFAULT: 0.
                type: int?
                inputBinding:
                  prefix: --shift
                  position: 1
              slocal:
                doc: |
                  SMALLLOCAL   The small nearby region in basepairs to calculate
                  dynamic lambda. This is used to capture the bias near
                  the peak summit region. Invalid if there is no control
                  data. If you set this to 0, MACS will skip slocal
                  lambda calculation. *Note* that MACS will always
                  perform a d-size local lambda calculation. The final
                  local bias should be the maximum of the lambda value
                  from d, slocal, and llocal size windows. DEFAULT: 1000
                type: int?
                inputBinding:
                  prefix: --slocal
                  position: 1
              to-large:
                doc: |-
                  When set, scale the small sample up to the bigger sample. By default, the bigger dataset will be scaled down towards the smaller dataset, which will lead to smaller p/qvalues and more specific results. Keep in mind that scaling down will bring down background noise more. DEFAULT: False 
                type: boolean?
                inputBinding:
                  prefix: --to-large
                  position: 1
              trackline:
                doc: |-
                  Tells MACS to include trackline with bedGraph files. To include this trackline while displaying bedGraph at UCSC genome browser, can show name and description of the file as well. However my suggestion is to convert bedGraph to bigWig, then show the smaller and faster binary bigWig file at UCSC genome browser, as well as downstream analysis. Require --bdg to be set. Default: Not include trackline. 
                type: boolean?
                inputBinding:
                  prefix: --trackline
                  position: 1
              treatment:
                doc: |-
                  Treatment sample file(s). If multiple files are given as -t A B C, then they will all be read and pooled together. IMPORTANT: the first sample will be used as the outputs basename.
                type: File[]
                inputBinding:
                  prefix: --treatment
                  position: 2
              verbose:
                doc: |
                  VERBOSE_LEVEL     Set verbose level of runtime message. 0: only show
                  critical message, 1: show additional warning message,
                  2: show process information, 3: show debug messages.
                  DEFAULT:2
                type: int?
                inputBinding:
                  prefix: --verbose
                  position: 1

            outputs:
              output_ext_frag_bdg_file:
                doc: Bedgraph with extended fragment pileup.
                type: File?
                outputBinding:
                  glob: |-
                    $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '_treat_pileup.bdg')
              output_peak_file:
                doc: Peak calling output file in narrowPeak|broadPeak format.
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '_peaks.*Peak')
                  outputEval: $(self[0])
              output_peak_summits_file:
                doc: Peaks summits bedfile.
                type: File?
                outputBinding:
                  glob: |-
                    $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '_summits.bed')
              output_peak_xls_file:
                doc: Peaks information/report file.
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '_peaks.xls')

            baseCommand:
            - macs2
            - callpeak
            arguments:
            - prefix: -n
              position: 1
              valueFrom: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
                ''))

            hints:
              DockerRequirement:
                dockerPull: dukegcb/macs2
          out:
          - output_peak_file
          - output_peak_summits_file
          - output_ext_frag_bdg_file
          - output_peak_xls_file
        peaks-bed-to-bigbed:
          in:
            type:
              valueFrom: bed6+4
            as: as_narrowPeak_file
            bed: trunk-peak-score/trunked_scores_peaks
            genome_sizes: input_genome_sizes
          scatter: bed
          run:
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
                doc: "-blockSize=N - Number of items to bundle in r-tree.  Default\
                  \ 256\n"
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
                doc: "-itemsPerSlot=N - Number of data points bundled at lowest level.\
                  \ Default 512\n"
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
          out:
          - bigbed
        spp:
          in:
            input_bam: input_bam_files
            nthreads: nthreads
            savp:
              valueFrom: ${return true}
          scatter:
          - input_bam
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
                doc: -fdr=<falseDisoveryRate> , false discovery rate threshold for
                  peak calling
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
                doc: 'overwrite (force remove) output files in case they exist. Default:
                  true'
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
                doc: -savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed
                  width peaks)
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
          out:
          - output_spp_cross_corr
          - output_spp_cross_corr_plot
        trunk-peak-score:
          in:
            peaks: peak-calling/output_peak_file
          scatter: peaks
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Trunk scores in ENCODE bed6+4 files

            inputs:
              peaks:
                type: File
                inputBinding:
                  position: 10000
              sep:
                type: string
                default: \t
                inputBinding:
                  prefix: -F
                  position: 2

            outputs:
              trunked_scores_peaks:
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.peaks.path.replace(/^.*[\\\/]/, '').replace(/\.([^/.]+)$/, "\.trunked_scores\.$1"))
            stdout: |-
              $(inputs.peaks.path.replace(/^.*[\\\/]/, '').replace(/\.([^/.]+)$/, "\.trunked_scores\.$1"))

            baseCommand: awk
            arguments:
            - position: 3
              valueFrom: BEGIN{OFS=FS}$5>1000{$5=1000}{print}

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - trunked_scores_peaks
    out:
    - output_spp_x_cross_corr
    - output_spp_cross_corr_plot
    - output_read_in_peak_count_within_replicate
    - output_peak_file
    - output_peak_bigbed_file
    - output_peak_summits_file
    - output_extended_peak_file
    - output_peak_xls_file
    - output_filtered_read_count_file
    - output_peak_count_within_replicate
  qc:
    in:
      default_adapters_file: default_adapters_file
      input_fastq_files: input_fastq_files
      nthreads: nthreads_qc
    run:
      cwlVersion: v1.0
      class: Workflow
      doc: 'ATAC-seq 01 QC - reads: SE'

      requirements:
      - class: ScatterFeatureRequirement
      - class: StepInputExpressionRequirement
      - class: InlineJavascriptRequirement

      inputs:
        default_adapters_file:
          doc: Adapters file
          type: File
        input_fastq_files:
          doc: Input fastq files
          type: File[]
        nthreads:
          doc: Number of threads.
          type: int

      outputs:
        output_count_raw_reads:
          type: File[]
          outputSource: count_raw_reads/output_read_count
        output_custom_adapters:
          type: File[]
          outputSource: overrepresented_sequence_extract/output_custom_adapters
        output_diff_counts:
          type: File[]
          outputSource: compare_read_counts/result
        output_fastqc_data_files:
          doc: FastQC data files
          type: File[]
          outputSource: extract_fastqc_data/output_fastqc_data_file
        output_fastqc_report_files:
          doc: FastQC reports in zip format
          type: File[]
          outputSource: fastqc/output_qc_report_file

      steps:
        compare_read_counts:
          in:
            file1: count_raw_reads/output_read_count
            file2: count_fastqc_reads/output_fastqc_read_count
          scatter:
          - file1
          - file2
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Compares 2 files

            inputs:
              brief:
                type: boolean
                default: true
                inputBinding:
                  prefix: --brief
                  position: 3
              file1:
                type: File
                inputBinding:
                  position: 1
              file2:
                type: File
                inputBinding:
                  position: 2

            outputs:
              result:
                type: File
                outputBinding:
                  glob: stdout.txt
            stdout: stdout.txt

            baseCommand: diff

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - result
        count_fastqc_reads:
          in:
            input_basename: extract_basename/output_basename
            input_fastqc_data: extract_fastqc_data/output_fastqc_data_file
          scatter:
          - input_fastqc_data
          - input_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts read count from fastqc_data.txt

            inputs:
              input_basename:
                type: string
              input_fastqc_data:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_fastqc_read_count:
                type: File
                outputBinding:
                  glob: $(inputs.input_basename + '.fastqc-read_count.txt')
            stdout: $(inputs.input_basename + '.fastqc-read_count.txt')

            baseCommand: count-fastqc_data-reads.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_fastqc_read_count
        count_raw_reads:
          in:
            input_basename: extract_basename/output_basename
            input_fastq_file: input_fastq_files
          scatter:
          - input_fastq_file
          - input_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Counts reads in a fastq file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_basename:
                type: string
              input_fastq_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_read_count:
                type: File
                outputBinding:
                  glob: $(inputs.input_basename + '.read_count.txt')
            stdout: $(inputs.input_basename + '.read_count.txt')

            baseCommand: count-fastq-reads.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_read_count
        extract_basename:
          in:
            input_file: input_fastq_files
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts the base name of a file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_basename:
                type: string
                outputBinding:
                  outputEval: |-
                    $(inputs.input_file.path.substr(inputs.input_file.path.lastIndexOf('/') + 1, inputs.input_file.path.lastIndexOf('.') - (inputs.input_file.path.lastIndexOf('/') + 1)))

            baseCommand: echo

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_basename
        extract_fastqc_data:
          in:
            input_basename: extract_basename/output_basename
            input_qc_report_file: fastqc/output_qc_report_file
          scatter:
          - input_qc_report_file
          - input_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |-
              Unzips a zipped fastqc report and returns the fastqc_data.txt file. Unzips the file to pipe and uses redirection

            inputs:
              extract_pattern:
                type: string
                default: '*/fastqc_data.txt'
                inputBinding:
                  position: 3
              input_basename:
                type: string
              input_qc_report_file:
                type: File
                inputBinding:
                  position: 2
              pipe:
                type: boolean
                default: true
                inputBinding:
                  prefix: -p
                  position: 1

            outputs:
              output_fastqc_data_file:
                type: File
                outputBinding:
                  glob: $(inputs.input_basename + '.fastqc_data.txt')
            stdout: $(inputs.input_basename + '.fastqc_data.txt')

            baseCommand: unzip

            hints:
              DockerRequirement:
                dockerPull: dukegcb/fastqc
          out:
          - output_fastqc_data_file
        fastqc:
          in:
            input_fastq_file: input_fastq_files
            threads: nthreads
          scatter: input_fastq_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool

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
          out:
          - output_qc_report_file
        overrepresented_sequence_extract:
          in:
            default_adapters_file: default_adapters_file
            input_basename: extract_basename/output_basename
            input_fastqc_data: extract_fastqc_data/output_fastqc_data_file
          scatter:
          - input_fastqc_data
          - input_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool

            inputs:
              default_adapters_file:
                doc: Adapters file in fasta format
                type: File
                inputBinding:
                  position: 2
              input_basename:
                doc: Name of the sample - used as a base name for generating output
                  files
                type: string
              input_fastqc_data:
                doc: fastqc_data.txt file from a fastqc report
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_custom_adapters:
                type: File
                outputBinding:
                  glob: $(inputs.input_basename + '.custom_adapters.fasta')

            baseCommand: overrepresented_sequence_extract.py
            arguments:
            - position: 3
              valueFrom: $(inputs.input_basename + '.custom_adapters.fasta')

            hints:
              DockerRequirement:
                dockerPull: reddylab/overrepresented_sequence_extract:1.0
          out:
          - output_custom_adapters
    out:
    - output_count_raw_reads
    - output_diff_counts
    - output_fastqc_report_files
    - output_fastqc_data_files
    - output_custom_adapters
  quant:
    in:
      input_bam_files: map/output_data_sorted_dedup_bam_files
      input_genome_sizes: genome_sizes_file
      nthreads: nthreads_quant
    run:
      cwlVersion: v1.0
      class: Workflow
      doc: ATAC-seq - Quantification

      requirements:
      - class: ScatterFeatureRequirement
      - class: StepInputExpressionRequirement
      - class: InlineJavascriptRequirement

      inputs:
        input_bam_files:
          type: File[]
        input_genome_sizes:
          type: File
        nthreads:
          type: int
          default: 1

      outputs:
        bigwig_norm_files:
          doc: signal files of pileup reads in RPKM
          type: File[]
          outputSource: bamcoverage/output_bam_coverage
        bigwig_raw_files:
          doc: Raw reads bigWig (signal) files
          type: File[]
          outputSource: bdg2bw-raw/output_bigwig

      steps:
        bamcoverage:
          in:
            bam: input_bam_files
            binSize:
              valueFrom: ${return 1}
            extendReads:
              valueFrom: ${return 200}
            normalizeUsing:
              valueFrom: RPKM
            numberOfProcessors: nthreads
            output_suffix:
              valueFrom: .rpkm.bw
          scatter: bam
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |
              usage: An example usage is:$ bamCoverage -b reads.bam -o coverage.bw

              This tool takes an alignment of reads or fragments as input (BAM file) and
              generates a coverage track (bigWig or bedGraph) as output. The coverage is
              calculated as the number of reads per bin, where bins are short consecutive
              counting windows of a defined size. It is possible to extended the length of
              the reads to better reflect the actual fragment length. *bamCoverage* offers
              normalization by scaling factor, Reads Per Kilobase per Million mapped reads
              (RPKM), and 1x depth (reads per genome coverage, RPGC).
              Required arguments:
                --bam BAM file, -b BAM file
                                      BAM file to process (default: None)
              Output:
                --outFileName FILENAME, -o FILENAME
                                      Output file name. (default: None)
                --outFileFormat {bigwig,bedgraph}, -of {bigwig,bedgraph}
                                      Output file type. Either "bigwig" or "bedgraph".
                                      (default: bigwig)
              Optional arguments:
                --help, -h            show this help message and exit
                --scaleFactor SCALEFACTOR
                                      The smooth length defines a window, larger than the
                                      binSize, to average the number of reads. For example,
                                      if the –binSize is set to 20 and the –smoothLength is
                                      set to 60, then, for each bin, the average of the bin
                                      and its left and right neighbors is considered.
                                      Any value smaller than –binSize will be ignored and
                                      no smoothing will be applied. (default: 1.0)
                --MNase               Determine nucleosome positions from MNase-seq data.
                                      Only 3 nucleotides at the center of each fragment are
                                      counted. The fragment ends are defined by the two mate
                                      reads. Only fragment lengthsbetween 130 - 200 bp are
                                      considered to avoid dinucleosomes or other
                                      artifacts.*NOTE*: Requires paired-end data. A bin size
                                      of 1 is recommended. (default: False)
                --filterRNAstrand {forward,reverse}
                                      Selects RNA-seq reads (single-end or paired-end) in
                                      the given strand. (default: None)
                --version             show program's version number and exit
                --binSize INT bp, -bs INT bp
                                      Size of the bins, in bases, for the output of the
                                      bigwig/bedgraph file. (default: 50)
                --region CHR:START:END, -r CHR:START:END
                                      Region of the genome to limit the operation to - this
                                      is useful when testing parameters to reduce the
                                      computing time. The format is chr:start:end, for
                                      example --region chr10 or --region
                                      chr10:456700:891000. (default: None)
                --blackListFileName BED file, -bl BED file
                                      A BED file containing regions that should be excluded
                                      from all analyses. Currently this works by rejecting
                                      genomic chunks that happen to overlap an entry.
                                      Consequently, for BAM files, if a read partially
                                      overlaps a blacklisted region or a fragment spans over
                                      it, then the read/fragment might still be considered.
                                      (default: None)
                --numberOfProcessors INT, -p INT
                                      Number of processors to use. Type "max/2" to use half
                                      the maximum number of processors or "max" to use all
                                      available processors. (default: max/2)
                --verbose, -v         Set to see processing messages. (default: False)
              Read coverage normalization options:
                --normalizeTo1x EFFECTIVE GENOME SIZE LENGTH
                                      Report read coverage normalized to 1x sequencing depth
                                      (also known as Reads Per Genomic Content (RPGC)).
                                      Sequencing depth is defined as: (total number of
                                      mapped reads * fragment length) / effective genome
                                      size. The scaling factor used is the inverse of the
                                      sequencing depth computed for the sample to match the
                                      1x coverage. To use this option, the effective genome
                                      size has to be indicated after the option. The
                                      effective genome size is the portion of the genome
                                      that is mappable. Large fractions of the genome are
                                      stretches of NNNN that should be discarded. Also, if
                                      repetitive regions were not included in the mapping of
                                      reads, the effective genome size needs to be adjusted
                                      accordingly. Common values are: mm9: 2,150,570,000;
                                      hg19:2,451,960,000; dm3:121,400,000 and
                                      ce10:93,260,000. See Table 2 of http://www.plosone.org
                                      /article/info:doi/10.1371/journal.pone.0030377 or http
                                      ://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_
                                      T1.html for several effective genome sizes. (default:
                                      None)
              --ignoreForNormalization IGNOREFORNORMALIZATION [IGNOREFORNORMALIZATION ...]
                                      A list of space-delimited chromosome names containing
                                      those chromosomes that should be excluded for
                                      computing the normalization. This is useful when
                                      considering samples with unequal coverage across
                                      chromosomes, like male samples. An usage examples is
                                      --ignoreForNormalization chrX chrM. (default: None)
                --skipNonCoveredRegions, --skipNAs
                                      This parameter determines if non-covered regions
                                      (regions without overlapping reads) in a BAM file
                                      should be skipped. The default is to treat those
                                      regions as having a value of zero. The decision to
                                      skip non-covered regions depends on the interpretation
                                      of the data. Non-covered regions may represent, for
                                      example, repetitive regions that should be skipped.
                                      (default: False)
                --smoothLength INT bp
                                      The smooth length defines a window, larger than the
                                      binSize, to average the number of reads. For example,
                                      if the --binSize is set to 20 and the --smoothLength
                                      is set to 60, then, for each bin, the average of the
                                      bin and its left and right neighbors is considered.
                                      Any value smaller than --binSize will be ignored and
                                      no smoothing will be applied. (default: None)
              Read processing options:
                --extendReads [INT bp], -e [INT bp]
                                      This parameter allows the extension of reads to
                                      fragment size. If set, each read is extended, without
                                      exception. *NOTE*: This feature is generally NOT
                                      recommended for spliced-read data, such as RNA-seq, as
                                      it would extend reads over skipped regions. *Single-
                                      end*: Requires a user specified value for the final
                                      fragment length. Reads that already exceed this
                                      fragment length will not be extended. *Paired-end*:
                                      Reads with mates are always extended to match the
                                      fragment size defined by the two read mates. Unmated
                                      reads, mate reads that map too far apart (>4x fragment
                                      length) or even map to different chromosomes are
                                      treated like single-end reads. The input of a fragment
                                      length value is optional. If no value is specified, it
                                      is estimated from the data (mean of the fragment size
                                      of all mate reads). (default: False)
                --ignoreDuplicates    If set, reads that have the same orientation and start
                                      position will be considered only once. If reads are
                                      paired, the mate's position also has to coincide to
                                      ignore a read. (default: False)
                --minMappingQuality INT
                                      If set, only reads that have a mapping quality score
                                      of at least this are considered. (default: None)
                --centerReads         By adding this option, reads are centered with respect
                                      to the fragment length. For paired-end data, the read
                                      is centered at the fragment length defined by the two
                                      ends of the fragment. For single-end data, the given
                                      fragment length is used. This option is useful to get
                                      a sharper signal around enriched regions. (default:
                                      False)
                --samFlagInclude INT  Include reads based on the SAM flag. For example, to
                                      get only reads that are the first mate, use a flag of
                                      64. This is useful to count properly paired reads only
                                      once, as otherwise the second mate will be also
                                      considered for the coverage. (default: None)
                --samFlagExclude INT  Exclude reads based on the SAM flag. For example, to
                                      get only reads that map to the forward strand, use
                                      --samFlagExclude 16, where 16 is the SAM flag for
                                      reads that map to the reverse strand. (default: None)

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              MNase:
                doc: |
                  Determine nucleosome positions from MNase-seq data.
                  Only 3 nucleotides at the center of each fragment are
                  counted. The fragment ends are defined by the two mate
                  reads. Only fragment lengthsbetween 130 - 200 bp are
                  considered to avoid dinucleosomes or other
                  artifacts.*NOTE*: Requires paired-end data. A bin size
                  of 1 is recommended. (default: False)
                type: boolean?
                inputBinding:
                  prefix: --MNase
                  position: 1
              bam:
                doc: 'BAM file to process '
                type: File
                secondaryFiles:
                - .bai
                inputBinding:
                  prefix: --bam
                  position: 1
              binSize:
                doc: |
                  INT bp
                  Size of the bins, in bases, for the output of the
                  bigwig/bedgraph file. (default: 50)
                type: int?
                inputBinding:
                  prefix: --binSize
                  position: 1
              blackListFileName:
                doc: |
                  BED file
                  A BED file containing regions that should be excluded
                  from all analyses. Currently this works by rejecting
                  genomic chunks that happen to overlap an entry.
                  Consequently, for BAM files, if a read partially
                  overlaps a blacklisted region or a fragment spans over
                  it, then the read/fragment might still be considered.
                  (default: None)
                type: File?
                inputBinding:
                  prefix: --blackListFileName
                  position: 1
              centerReads:
                doc: |
                  By adding this option, reads are centered with respect
                  to the fragment length. For paired-end data, the read
                  is centered at the fragment length defined by the two
                  ends of the fragment. For single-end data, the given
                  fragment length is used. This option is useful to get
                  a sharper signal around enriched regions. (default:
                  False)
                type: boolean?
                inputBinding:
                  prefix: --centerReads
                  position: 1
              extendReads:
                doc: |
                  INT bp
                  This parameter allows the extension of reads to
                  fragment size. If set, each read is extended, without
                  exception. *NOTE*: This feature is generally NOT
                  recommended for spliced-read data, such as RNA-seq, as
                  it would extend reads over skipped regions. *Single-
                  end*: Requires a user specified value for the final
                  fragment length. Reads that already exceed this
                  fragment length will not be extended. *Paired-end*:
                  Reads with mates are always extended to match the
                  fragment size defined by the two read mates. Unmated
                  reads, mate reads that map too far apart (>4x fragment
                  length) or even map to different chromosomes are
                  treated like single-end reads. The input of a fragment
                  length value is optional. If no value is specified, it
                  is estimated from the data (mean of the fragment size
                  of all mate reads). (default: False)
                type: int?
                inputBinding:
                  prefix: --extendReads
                  position: 1
              filterRNAstrand:
                doc: |
                  {forward,reverse}
                  Selects RNA-seq reads (single-end or paired-end) in
                  the given strand. (default: None)
                type: string?
                inputBinding:
                  prefix: --filterRNAstrand
                  position: 1
              ignoreDuplicates:
                doc: |
                  If set, reads that have the same orientation and start
                  position will be considered only once. If reads are
                  paired, the mate's position also has to coincide to
                  ignore a read. (default: False)
                type: boolean?
                inputBinding:
                  prefix: --ignoreDuplicates
                  position: 1
              ignoreForNormalization:
                doc: |
                  A list of space-delimited chromosome names containing those chromosomes
                  that should be excluded for computing the normalization. This is useful
                  when considering samples with unequal coverage across chromosomes, like
                  male samples. An usage examples is –ignoreForNormalization chrX chrM.
                type:
                - 'null'
                - type: array
                  items: string
                inputBinding:
                  prefix: --ignoreForNormalization
                  position: 1
              minMappingQuality:
                doc: |
                  INT
                  If set, only reads that have a mapping quality score
                  of at least this are considered. (default: None)
                type: int?
                inputBinding:
                  prefix: --minMappingQuality
                  position: 1
              normalizeUsing:
                doc: |
                  Possible choices: RPKM, CPM, BPM, RPGC

                  Use one of the entered methods to normalize the number of reads per bin.
                  By default, no normalization is performed. RPKM = Reads Per Kilobase per
                  Million mapped reads; CPM = Counts Per Million mapped reads, same as
                  CPM in RNA-seq; BPM = Bins Per Million mapped reads, same as TPM in
                  RNA-seq; RPGC = reads per genomic content (1x normalization);
                  Mapped reads are considered after blacklist filtering (if applied).
                  RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)).
                  CPM (per bin) = number of reads per bin / number of mapped reads (in millions).
                  BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions).
                  RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage.
                  This scaling factor, in turn, is determined from the sequencing depth:
                  (total number of mapped reads * fragment length) / effective genome size.
                  The scaling factor used is the inverse of the sequencing depth
                  computed for the sample to match the 1x coverage.
                  This option requires –effectiveGenomeSize.
                  Each read is considered independently, if you want to only count one
                  mate from a pair in paired-end data, then use
                  the –samFlagInclude/–samFlagExclude options.
                type: string?
                inputBinding:
                  prefix: --normalizeUsing
                  position: 1
              numberOfProcessors:
                doc: |
                  INT
                  Number of processors to use. Type "max/2" to use half
                  the maximum number of processors or "max" to use all
                  available processors. (default: max/2)
                type: int?
                inputBinding:
                  prefix: --numberOfProcessors
                  position: 1
              outFileFormat:
                doc: |
                  {bigwig,bedgraph}, -of {bigwig,bedgraph}
                  Output file type. Either "bigwig" or "bedgraph".
                  (default: bigwig)
                type: string
                default: bigwig
                inputBinding:
                  prefix: --outFileFormat
                  position: 1
              outFileName:
                doc: |
                  FILENAME
                  Output file name. (default: input BAM filename with bigwig [*.bw] or bedgraph [*.bdg] extension.)
                type: string?
              output_suffix:
                doc: Suffix used for output file (input BAM filename + suffix)
                type: string?
              region:
                doc: |
                  CHR:START:END
                  Region of the genome to limit the operation to - this
                  is useful when testing parameters to reduce the
                  computing time. The format is chr:start:end, for
                  example --region chr10 or --region
                  chr10:456700:891000. (default: None)
                type: string?
                inputBinding:
                  prefix: --region
                  position: 1
              samFlagExclude:
                doc: |
                  INT  
                  Exclude reads based on the SAM flag. For example, to
                  get only reads that map to the forward strand, use
                  --samFlagExclude 16, where 16 is the SAM flag for
                  reads that map to the reverse strand. (default: None)
                type: int?
                inputBinding:
                  prefix: --samFlagExclude
                  position: 1
              samFlagInclude:
                doc: |
                  INT  
                  Include reads based on the SAM flag. For example, to
                  get only reads that are the first mate, use a flag of
                  64. This is useful to count properly paired reads only
                  once, as otherwise the second mate will be also
                  considered for the coverage. (default: None)
                type: int?
                inputBinding:
                  prefix: --samFlagInclude
                  position: 1
              scaleFactor:
                doc: |
                  The smooth length defines a window, larger than the
                  binSize, to average the number of reads. For example,
                  if the –binSize is set to 20 and the –smoothLength is
                  set to 60, then, for each bin, the average of the bin
                  and its left and right neighbors is considered.
                  Any value smaller than –binSize will be ignored and
                  no smoothing will be applied. (default: 1.0)
                type: float?
                inputBinding:
                  prefix: --scaleFactor
                  position: 1
              skipNonCoveredRegions:
                doc: |
                  --skipNonCoveredRegions, --skipNAs
                  This parameter determines if non-covered regions
                  (regions without overlapping reads) in a BAM file
                  should be skipped. The default is to treat those
                  regions as having a value of zero. The decision to
                  skip non-covered regions depends on the interpretation
                  of the data. Non-covered regions may represent, for
                  example, repetitive regions that should be skipped.
                  (default: False)
                type: boolean?
                inputBinding:
                  prefix: --skipNonCoveredRegions
                  position: 1
              smoothLength:
                doc: |
                  INT bp
                  The smooth length defines a window, larger than the
                  binSize, to average the number of reads. For example,
                  if the --binSize is set to 20 and the --smoothLength
                  is set to 60, then, for each bin, the average of the
                  bin and its left and right neighbors is considered.
                  Any value smaller than --binSize will be ignored and
                  no smoothing will be applied. (default: None)
                  Read processing options:
                type: int?
                inputBinding:
                  prefix: --smoothLength
                  position: 1
              verbose:
                doc: "--verbose         \nSet to see processing messages. (default:\
                  \ False)\n"
                type: boolean?
                inputBinding:
                  prefix: --verbose
                  position: 1
              version:
                doc: show program's version number and exit
                type: boolean?
                inputBinding:
                  prefix: --version
                  position: 1

            outputs:
              output_bam_coverage:
                type: File
                outputBinding:
                  glob: |-
                    ${ if (inputs.outFileName) return inputs.outFileName; if (inputs.output_suffix) return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + inputs.output_suffix; if (inputs.outFileFormat == "bedgraph") return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".bdg"; return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".bw"; }

            baseCommand: bamCoverage
            arguments:
            - prefix: --outFileName
              position: 3
              valueFrom: |-
                ${ if (inputs.outFileName) return inputs.outFileName; if (inputs.output_suffix) return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + inputs.output_suffix; if (inputs.outFileFormat == "bedgraph") return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".bdg"; return inputs.bam.path.replace(/^.*[\\\/]/, "").replace(/\.[^/.]+$/, "") + ".bw"; }

            hints:
              DockerRequirement:
                dockerPull: reddylab/deeptools:3.0.1
          out:
          - output_bam_coverage
        bdg2bw-raw:
          in:
            bed_graph: bedsort_genomecov/bed_file_sorted
            genome_sizes: input_genome_sizes
            output_suffix:
              valueFrom: .raw.bw
          scatter: bed_graph
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: 'Tool:   bedGraphToBigWig v 4 - Convert a bedGraph file to bigWig
              format.'

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              bed_graph:
                doc: "\tbed_graph is a four column file in the format: <chrom> <start>\
                  \ <end> <value>\n"
                type: File
                inputBinding:
                  position: 1
              genome_sizes:
                doc: "\tgenome_sizes is two column: <chromosome name> <size in bases>.\n"
                type: File
                inputBinding:
                  position: 2
              output_suffix:
                type: string
                default: .bw

            outputs:
              output_bigwig:
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.bed_graph.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

            baseCommand: bedGraphToBigWig
            arguments:
            - position: 3
              valueFrom: |-
                $(inputs.bed_graph.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, "") + inputs.output_suffix)

            hints:
              DockerRequirement:
                dockerPull: dukegcb/bedgraphtobigwig
          out:
          - output_bigwig
        bedsort_genomecov:
          in:
            bed_file: bedtools_genomecov/output_bedfile
          scatter: bed_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |
              bedSort - Sort a .bed file by chrom,chromStart
              usage:
                 bedSort in.bed out.bed
              in.bed and out.bed may be the same.

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              bed_file:
                doc: Bed or bedGraph file to be sorted
                type: File
                inputBinding:
                  position: 1

            outputs:
              bed_file_sorted:
                type: File
                outputBinding:
                  glob: $(inputs.bed_file.path.replace(/^.*[\\\/]/, '') + "_sorted")

            baseCommand: bedSort
            arguments:
            - position: 2
              valueFrom: $(inputs.bed_file.path.replace(/^.*[\\\/]/, '') + "_sorted")

            hints:
              DockerRequirement:
                dockerPull: dleehr/docker-hubutils
          out:
          - bed_file_sorted
        bedtools_genomecov:
          in:
            bg:
              valueFrom: ${return true}
            g: input_genome_sizes
            ibam: input_bam_files
          scatter: ibam
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |-
              Tool:    bedtools genomecov (aka genomeCoverageBed)
              Version: v2.25.0
              Summary: Compute the coverage of a feature file among a genome.

              Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>

              Options: 
              	-ibam		The input file is in BAM format.
              			Note: BAM _must_ be sorted by position

              	-d		Report the depth at each genome position (with one-based coordinates).
              			Default behavior is to report a histogram.

              	-dz		Report the depth at each genome position (with zero-based coordinates).
              			Reports only non-zero positions.
              			Default behavior is to report a histogram.

              	-bg		Report depth in BedGraph format. For details, see:
              			genome.ucsc.edu/goldenPath/help/bedgraph.html

              	-bga		Report depth in BedGraph format, as above (-bg).
              			However with this option, regions with zero 
              			coverage are also reported. This allows one to
              			quickly extract all regions of a genome with 0 
              			coverage by applying: "grep -w 0$" to the output.

              	-split		Treat "split" BAM or BED12 entries as distinct BED intervals.
              			when computing coverage.
              			For BAM files, this uses the CIGAR "N" and "D" operations 
              			to infer the blocks for computing coverage.
              			For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
              			fields (i.e., columns 10,11,12).

              	-strand		Calculate coverage of intervals from a specific strand.
              			With BED files, requires at least 6 columns (strand is column 6). 
              			- (STRING): can be + or -

              	-5		Calculate coverage of 5" positions (instead of entire interval).

              	-3		Calculate coverage of 3" positions (instead of entire interval).

              	-max		Combine all positions with a depth >= max into
              			a single bin in the histogram. Irrelevant
              			for -d and -bedGraph
              			- (INTEGER)

              	-scale		Scale the coverage by a constant factor.
              			Each coverage value is multiplied by this factor before being reported.
              			Useful for normalizing coverage by, e.g., reads per million (RPM).
              			- Default is 1.0; i.e., unscaled.
              			- (FLOAT)

              	-trackline	Adds a UCSC/Genome-Browser track line definition in the first line of the output.
              			- See here for more details about track line definition:
              			      http://genome.ucsc.edu/goldenPath/help/bedgraph.html
              			- NOTE: When adding a trackline definition, the output BedGraph can be easily
              			      uploaded to the Genome Browser as a custom track,
              			      BUT CAN NOT be converted into a BigWig file (w/o removing the first line).

              	-trackopts	Writes additional track line definition parameters in the first line.
              			- Example:
              			   -trackopts 'name="My Track" visibility=2 color=255,30,30'
              			   Note the use of single-quotes if you have spaces in your parameters.
              			- (TEXT)

              Notes: 
              	(1) The genome file should tab delimited and structured as follows:
              	 <chromName><TAB><chromSize>

              	For example, Human (hg19):
              	chr1	249250621
              	chr2	243199373
              	...
              	chr18_gl000207_random	4262

              	(2) The input BED (-i) file must be grouped by chromosome.
              	 A simple "sort -k 1,1 <BED> > <BED>.sorted" will suffice.

              	(3) The input BAM (-ibam) file must be sorted by position.
              	 A "samtools sort <BAM>" should suffice.

              Tips: 
              	One can use the UCSC Genome Browser's MySQL database to extract
              	chromosome sizes. For example, H. sapiens:

              	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
              	"select chrom, size from hg19.chromInfo" > hg19.genome

            requirements:
              InlineJavascriptRequirement: {}
              ShellCommandRequirement: {}

            inputs:
              '3':
                doc: "\tCalculate coverage of 3\" positions (instead of entire interval).\n"
                type: boolean?
                inputBinding:
                  prefix: '-3'
                  position: 1
              '5':
                doc: "\tCalculate coverage of 5\" positions (instead of entire interval).\n"
                type: boolean?
                inputBinding:
                  prefix: '-5'
                  position: 1
              bg:
                doc: |
                  	Report depth in BedGraph format. For details, see:
                  genome.ucsc.edu/goldenPath/help/bedgraph.html
                type: boolean?
                inputBinding:
                  prefix: -bg
                  position: 1
              bga:
                doc: |
                  	Report depth in BedGraph format, as above (-bg).
                  However with this option, regions with zero
                  coverage are also reported. This allows one to
                  quickly extract all regions of a genome with 0
                  coverage by applying: "grep -w 0$" to the output.
                type: boolean?
                inputBinding:
                  prefix: -bga
                  position: 1
              d:
                doc: |
                  	Report the depth at each genome position (with one-based coordinates).
                  Default behavior is to report a histogram.
                type: boolean?
                inputBinding:
                  prefix: -d
                  position: 1
              dz:
                doc: |
                  	Report the depth at each genome position (with zero-based coordinates).
                  Reports only non-zero positions.
                  Default behavior is to report a histogram.
                type: boolean?
                inputBinding:
                  prefix: -dz
                  position: 1
              g:
                doc: <genome sizes>
                type: File
                inputBinding:
                  prefix: -g
                  position: 3
              ibam:
                doc: "\tThe input file is in BAM format.\nNote: BAM _must_ be sorted\
                  \ by position\n"
                type: File
                inputBinding:
                  prefix: -ibam
                  position: 2
              max:
                doc: |
                  	Combine all positions with a depth >= max into
                  a single bin in the histogram. Irrelevant
                  for -d and -bedGraph
                  - (INTEGER)
                type: int?
                inputBinding:
                  prefix: -max
                  position: 1
              scale:
                doc: |
                  	Scale the coverage by a constant factor.
                  Each coverage value is multiplied by this factor before being reported.
                  Useful for normalizing coverage by, e.g., reads per million (RPM).
                  - Default is 1.0; i.e., unscaled.
                  - (FLOAT)
                type: float?
                inputBinding:
                  prefix: -scale
                  position: 1
              split:
                doc: |
                  	Treat "split" BAM or BED12 entries as distinct BED intervals.
                  when computing coverage.
                  For BAM files, this uses the CIGAR "N" and "D" operations
                  to infer the blocks for computing coverage.
                  For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
                  fields (i.e., columns 10,11,12).
                type: boolean?
                inputBinding:
                  prefix: -split
                  position: 1
              strand:
                doc: |
                  	Calculate coverage of intervals from a specific strand.
                  With BED files, requires at least 6 columns (strand is column 6).
                  - (STRING): can be + or -
                type: string?
                inputBinding:
                  prefix: -strand
                  position: 1
              trackline:
                doc: |
                  Adds a UCSC/Genome-Browser track line definition in the first line of the output.
                  - See here for more details about track line definition:
                  http://genome.ucsc.edu/goldenPath/help/bedgraph.html
                  - NOTE: When adding a trackline definition, the output BedGraph can be easily
                  uploaded to the Genome Browser as a custom track,
                  BUT CAN NOT be converted into a BigWig file (w/o removing the first line).
                type: boolean?
                inputBinding:
                  prefix: -trackline
                  position: 1
              trackopts:
                doc: |
                  Writes additional track line definition parameters in the first line.
                  - Example:
                  -trackopts 'name="My Track" visibility=2 color=255,30,30'
                  Note the use of single-quotes if you have spaces in your parameters.
                  - (TEXT)
                type: string?
                inputBinding:
                  prefix: -trackopts
                  position: 1

            outputs:
              output_bedfile:
                type: File
                outputBinding:
                  glob: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
                    '') + '.bdg')
            stdout: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
              '') + '.bdg')

            baseCommand:
            - bedtools
            - genomecov

            hints:
              DockerRequirement:
                dockerPull: dukegcb/bedtools
          out:
          - output_bedfile
    out:
    - bigwig_raw_files
    - bigwig_norm_files
  trimm:
    in:
      input_adapters_files: qc/output_custom_adapters
      input_fastq_files: input_fastq_files
      nthreads: nthreads_trimm
      trimmomatic_jar_path: trimmomatic_jar_path
      trimmomatic_java_opts: trimmomatic_java_opts
    run:
      cwlVersion: v1.0
      class: Workflow
      doc: 'ATAC-seq 02 trimming - reads: SE'

      requirements:
      - class: ScatterFeatureRequirement
      - class: StepInputExpressionRequirement
      - class: InlineJavascriptRequirement

      inputs:
        input_adapters_files:
          doc: Input adapters files
          type: File[]
        input_fastq_files:
          doc: Input fastq files
          type: File[]
        nthreads:
          doc: Number of threads
          type: int
          default: 1
        quality_score:
          type: string
          default: -phred33
        trimmomatic_jar_path:
          doc: Trimmomatic Java jar file
          type: string
          default: /usr/share/java/trimmomatic.jar
        trimmomatic_java_opts:
          doc: JVM arguments should be a quoted, space separated list
          type: string?

      outputs:
        output_data_fastq_trimmed_files:
          doc: Trimmed fastq files
          type: File[]
          outputSource: trimmomatic/output_read1_trimmed_file
        output_trimmed_fastq_read_count:
          doc: Trimmed read counts of fastq files
          type: File[]
          outputSource: count_fastq_reads/output_read_count

      steps:
        count_fastq_reads:
          in:
            input_basename: extract_basename/output_basename
            input_fastq_file: trimmomatic/output_read1_trimmed_file
          scatter:
          - input_fastq_file
          - input_basename
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Counts reads in a fastq file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_basename:
                type: string
              input_fastq_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_read_count:
                type: File
                outputBinding:
                  glob: $(inputs.input_basename + '.read_count.txt')
            stdout: $(inputs.input_basename + '.read_count.txt')

            baseCommand: count-fastq-reads.sh

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_read_count
        extract_basename:
          in:
            input_file: trimmomatic/output_read1_trimmed_file
          scatter: input_file
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: Extracts the base name of a file

            requirements:
              InlineJavascriptRequirement: {}

            inputs:
              input_file:
                type: File
                inputBinding:
                  position: 1

            outputs:
              output_basename:
                type: string
                outputBinding:
                  outputEval: |-
                    $(inputs.input_file.path.substr(inputs.input_file.path.lastIndexOf('/') + 1, inputs.input_file.path.lastIndexOf('.') - (inputs.input_file.path.lastIndexOf('/') + 1)))

            baseCommand: echo

            hints:
              DockerRequirement:
                dockerPull: reddylab/workflow-utils:ggr
          out:
          - output_basename
        trimmomatic:
          in:
            end_mode:
              valueFrom: SE
            illuminaclip:
              valueFrom: 2:30:15
            input_adapters_file: input_adapters_files
            input_read1_fastq_file: input_fastq_files
            java_opts: trimmomatic_java_opts
            leading:
              valueFrom: ${return 3}
            minlen:
              valueFrom: ${return 15}
            nthreads: nthreads
            phred:
              valueFrom: '33'
            slidingwindow:
              valueFrom: 4:20
            trailing:
              valueFrom: ${return 3}
            trimmomatic_jar_path: trimmomatic_jar_path
          scatter:
          - input_read1_fastq_file
          - input_adapters_file
          scatterMethod: dotproduct
          run:
            cwlVersion: v1.0
            class: CommandLineTool
            doc: |
              Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
              Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
              depending on the library preparation and downstream application.
              There are two major modes of the program: Paired end mode and Single end mode. The
              paired end mode will maintain correspondence of read pairs and also use the additional
              information contained in paired reads to better find adapter or PCR primer fragments
              introduced by the library preparation process.
              Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
              depending on the Illumina pipeline used).

            requirements:
              InlineJavascriptRequirement: {}
              ShellCommandRequirement: {}

            inputs:
              avgqual:
                doc: |
                  <quality>
                  Drop the read if the average quality is below the specified level
                  <quality>: Specifies the minimum average quality required to keep a read.
                type: int?
                inputBinding:
                  prefix: 'AVGQUAL:'
                  position: 101
                  separate: false
              crop:
                doc: |
                  <length>
                  Removes bases regardless of quality from the end of the read, so that the read has maximally
                  the specified length after this step has been performed. Steps performed after CROP might of
                  course further shorten the read.
                  <length>: The number of bases to keep, from the start of the read.
                type: int?
                inputBinding:
                  prefix: 'CROP:'
                  position: 13
                  separate: false
              end_mode:
                doc: "SE|PE\nSingle End (SE) or Paired End (PE) mode\n"
                type: string
                inputBinding:
                  position: 3
              headcrop:
                doc: |
                  <length>
                  Removes the specified number of bases, regardless of quality, from the beginning of the read.
                  <length>: The number of bases to keep, from the start of the read.
                type: int?
                inputBinding:
                  prefix: 'HEADCROP:'
                  position: 13
                  separate: false
              illuminaclip:
                doc: |
                  <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
                  Find and remove Illumina adapters.
                  REQUIRED:
                  <fastaWithAdaptersEtc>: specifies the path to a fasta file containing all the adapters, PCR sequences etc.
                  The naming of the various sequences within this file determines how they are used. See below.
                  <seedMismatches>: specifies the maximum mismatch count which will still allow a full match to be performed
                  <palindromeClipThreshold>: specifies how accurate the match between the two 'adapter ligated' reads must be
                  for PE palindrome read alignment.
                  <simpleClipThreshold>: specifies how accurate the match between any adapter etc. sequence must be against a read
                  OPTIONAL:
                  <minAdapterLength>: In addition to the alignment score, palindrome mode can verify
                  that a minimum length of adapter has been detected. If unspecified, this defaults to 8 bases,
                  for historical reasons. However, since palindrome mode has a very low false positive rate, this
                  can be safely reduced, even down to 1, to allow shorter adapter fragments to be removed.
                  <keepBothReads>: After read-though has been detected by palindrome mode, and the
                  adapter sequence removed, the reverse read contains the same sequence information as the
                  forward read, albeit in reverse complement. For this reason, the default behaviour is to
                  entirely drop the reverse read. By specifying "true" for this parameter, the reverse read will
                  also be retained, which may be useful e.g. if the downstream tools cannot handle a
                  combination of paired and unpaired reads.
                type: string
              input_adapters_file:
                doc: |-
                  FASTA file containing adapters, PCR sequences, etc. It is used to search for and remove these sequences in the input FASTQ file(s)
                type: File
              input_read1_fastq_file:
                doc: FASTQ file for input read (read R1 in Paired End mode)
                type: File
                inputBinding:
                  position: 5
              input_read2_fastq_file:
                doc: FASTQ file for read R2 in Paired End mode
                type: File?
                inputBinding:
                  position: 6
              java_opts:
                doc: |-
                  JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")
                type: string?
                inputBinding:
                  position: 1
                  shellQuote: false
              leading:
                doc: |
                  <quality>
                  Remove low quality bases from the beginning. As long as a base has a value below this
                  threshold the base is removed and the next base will be investigated.
                  <quality>: Specifies the minimum quality required to keep a base.
                type: int?
                inputBinding:
                  prefix: 'LEADING:'
                  position: 14
                  separate: false
              log_filename:
                doc: |
                  <ouptut log file name>
                  Specifying a trimlog file creates a log of all read trimmings, indicating the following details:
                    the read name
                    the surviving sequence length
                    the location of the first surviving base, aka. the amount trimmed from the start
                    the location of the last surviving base in the original read
                    the amount trimmed from the end
                  <ouptut log file name>: filename for the generated output log file.
                type: string?
                inputBinding:
                  prefix: -trimlog
                  position: 4
              maxinfo:
                doc: |
                  <targetLength>:<strictness>
                  Performs an adaptive quality trim, balancing the benefits of retaining longer reads against the
                  costs of retaining bases with errors.
                  <targetLength>: This specifies the read length which is likely to allow the
                  location of the read within the target sequence to be determined.
                  <strictness>: This value, which should be set between 0 and 1, specifies the
                  balance between preserving as much read length as possible vs. removal of incorrect
                  bases. A low value of this parameter (<0.2) favours longer reads, while a high value
                  (>0.8) favours read correctness.
                type: string?
                inputBinding:
                  prefix: 'MAXINFO:'
                  position: 15
                  separate: false
              minlen:
                doc: |
                  <length>
                  This module removes reads that fall below the specified minimal length. If required, it should
                  normally be after all other processing steps. Reads removed by this step will be counted and
                  included in the "dropped reads" count presented in the trimmomatic summary.
                  <length>:  Specifies the minimum length of reads to be kept
                type: int?
                inputBinding:
                  prefix: 'MINLEN:'
                  position: 100
                  separate: false
              nthreads:
                doc: Number of threads
                type: int
                default: 1
                inputBinding:
                  prefix: -threads
                  position: 4
              phred:
                doc: |
                  "33"|"64"
                  -phred33 ("33") or -phred64 ("64") specifies the base quality encoding. Default: -phred64
                type: string
                default: '64'
                inputBinding:
                  prefix: -phred
                  position: 4
                  separate: false
              slidingwindow:
                doc: |
                  <windowSize>:<requiredQuality>
                  Perform a sliding window trimming, cutting once the average quality within the window falls
                  below a threshold. By considering multiple bases, a single poor quality base will not cause the
                  removal of high quality data later in the read.
                  <windowSize>: specifies the number of bases to average across
                  <requiredQuality>: specifies the average quality required
                type: string?
                inputBinding:
                  prefix: 'SLIDINGWINDOW:'
                  position: 15
                  separate: false
              tophred33:
                doc: This (re)encodes the quality part of the FASTQ file to base 33.
                type: boolean?
                inputBinding:
                  prefix: TOPHRED33
                  position: 12
                  separate: false
              tophred64:
                doc: This (re)encodes the quality part of the FASTQ file to base 64.
                type: boolean?
                inputBinding:
                  prefix: TOPHRED64
                  position: 12
                  separate: false
              trailing:
                doc: |
                  <quality>
                  Remove low quality bases from the end. As long as a base has a value below this threshold
                  the base is removed and the next base (which as trimmomatic is starting from the 3" prime end
                  would be base preceding the just removed base) will be investigated. This approach can be
                  used removing the special illumina "low quality segment" regions (which are marked with
                  quality score of 2), but we recommend Sliding Window or MaxInfo instead
                  <quality>: Specifies the minimum quality required to keep a base.
                type: int?
                inputBinding:
                  prefix: 'TRAILING:'
                  position: 14
                  separate: false
              trimmomatic_jar_path:
                type: string
                inputBinding:
                  prefix: -jar
                  position: 2

            outputs:
              output_log_file:
                doc: Trimmomatic Log file.
                type: File?
                outputBinding:
                  glob: $(inputs.log_filename)
              output_read1_trimmed_file:
                type: File
                outputBinding:
                  glob: |-
                    $(inputs.input_read1_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.fastq')
              output_read1_trimmed_unpaired_file:
                type: File?
                outputBinding:
                  glob: |
                    ${
                      if (inputs.end_mode == "PE")
                        return inputs.input_read1_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.unpaired.trimmed.fastq';
                      return null;
                    }
              output_read2_trimmed_paired_file:
                type: File?
                outputBinding:
                  glob: |
                    ${
                      if (inputs.end_mode == "PE" && inputs.input_read2_fastq_file)
                        return inputs.input_read2_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.fastq';
                      return null;
                    }
              output_read2_trimmed_unpaired_file:
                type: File?
                outputBinding:
                  glob: |
                    ${
                      if (inputs.end_mode == "PE" && inputs.input_read2_fastq_file)
                        return inputs.input_read2_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.unpaired.trimmed.fastq';
                      return null;
                    }

            baseCommand: java
            arguments:
            - position: 1
              valueFrom: $("-Djava.io.tmpdir="+runtime.tmpdir)
              shellQuote: false
            - position: 7
              valueFrom: |-
                $(inputs.input_read1_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.fastq')
            - position: 8
              valueFrom: |
                ${
                  if (inputs.end_mode == "PE" && inputs.input_read2_fastq_file)
                    return inputs.input_read1_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.unpaired.fastq';
                  return null;
                }
            - position: 9
              valueFrom: |
                ${
                  if (inputs.end_mode == "PE" && inputs.input_read2_fastq_file)
                    return inputs.input_read2_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.fastq';
                  return null;
                }
            - position: 10
              valueFrom: |
                ${
                  if (inputs.end_mode == "PE" && inputs.input_read2_fastq_file)
                    return inputs.input_read2_fastq_file.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.trimmed.unpaired.fastq';
                  return null;
                }
            - position: 11
              valueFrom: $("ILLUMINACLIP:" + inputs.input_adapters_file.path + ":"+
                inputs.illuminaclip)

            hints:
              DockerRequirement:
                dockerPull: dukegcb/trimmomatic
          out:
          - output_read1_trimmed_file
    out:
    - output_data_fastq_trimmed_files
    - output_trimmed_fastq_read_count
id: |-
  https://api.sbgenomics.com/v2/apps/kghosesbg/sbpla-31744/ATAC-seq-pipeline-se/2/raw/
sbg:appVersion:
- v1.0
sbg:content_hash: ad9474546d1d7aba5aa20e3c7a03b5429e5f8ec1d18be92cbab7315600a6bce48
sbg:contributors:
- kghosesbg
sbg:createdBy: kghosesbg
sbg:createdOn: 1580500895
sbg:id: kghosesbg/sbpla-31744/ATAC-seq-pipeline-se/2
sbg:image_url: |-
  https://igor.sbgenomics.com/ns/brood/images/kghosesbg/sbpla-31744/ATAC-seq-pipeline-se/2.png
sbg:latestRevision: 2
sbg:modifiedBy: kghosesbg
sbg:modifiedOn: 1581699121
sbg:project: kghosesbg/sbpla-31744
sbg:projectName: SBPLA-31744
sbg:publisher: sbg
sbg:revision: 2
sbg:revisionNotes: |-
  Uploaded using sbpack v2020.02.14. 
  Source: https://raw.githubusercontent.com/Duke-GCB/GGR-cwl/master/v1.0/ATAC-seq_pipeline/pipeline-se.cwl
sbg:revisionsInfo:
- sbg:modifiedBy: kghosesbg
  sbg:modifiedOn: 1580500895
  sbg:revision: 0
  sbg:revisionNotes: |-
    Uploaded using sbpack. Source: https://raw.githubusercontent.com/Duke-GCB/GGR-cwl/master/v1.0/ATAC-seq_pipeline/pipeline-se.cwl
- sbg:modifiedBy: kghosesbg
  sbg:modifiedOn: 1580742764
  sbg:revision: 1
  sbg:revisionNotes: Just moved a node
- sbg:modifiedBy: kghosesbg
  sbg:modifiedOn: 1581699121
  sbg:revision: 2
  sbg:revisionNotes: |-
    Uploaded using sbpack v2020.02.14. 
    Source: https://raw.githubusercontent.com/Duke-GCB/GGR-cwl/master/v1.0/ATAC-seq_pipeline/pipeline-se.cwl
sbg:sbgMaintained: false
sbg:validationErrors:
- 'Required input is not set: #qc.input_fastq_files'
- 'Required input is not set: #qc.default_adapters_file'
- 'Required input is not set: #qc.nthreads'
- 'Required input is not set: #trimm.input_fastq_files'
- 'Required input is not set: #trimm.input_adapters_files'
- 'Required input is not set: #map.input_fastq_files'
- 'Required input is not set: #map.genome_sizes_file'
- 'Required input is not set: #map.genome_ref_first_index_file'
- 'Required input is not set: #peak_call.input_bam_files'
- 'Required input is not set: #peak_call.input_genome_sizes'
- 'Required input is not set: #peak_call.as_narrowPeak_file'
- 'Required input is not set: #quant.input_bam_files'
- 'Required input is not set: #quant.input_genome_sizes'
