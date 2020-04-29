class: CommandLineTool
cwlVersion: v1.0
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
    doc: "\tThe input file is in BAM format.\nNote: BAM _must_ be sorted by position\n"
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
      glob: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') +
        '.bdg')
stdout: $(inputs.ibam.path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/, '') + '.bdg')

baseCommand:
- bedtools
- genomecov

hints:
  DockerRequirement:
    dockerPull: dukegcb/bedtools
