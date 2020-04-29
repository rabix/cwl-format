class: CommandLineTool
cwlVersion: v1.0
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
  valueFrom: $("ILLUMINACLIP:" + inputs.input_adapters_file.path + ":"+ inputs.illuminaclip)

hints:
  DockerRequirement:
    dockerPull: dukegcb/trimmomatic
