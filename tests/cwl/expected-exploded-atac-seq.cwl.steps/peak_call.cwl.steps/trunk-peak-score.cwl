#!/usr/bin/env cwl-runner

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
