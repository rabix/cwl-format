# CWL Format

CWL Format is a very opinionated code formatter for CWL.
It outputs CWL in a standardized YAML format. It has no settings or
options because you have better things to do with your time. And because
cwl-format is always correct.

This repository lists the formatting rules and also contains a Python
implementation of the formatter.

## Rules

- Comments are not preserved **Do not use this if comments in the YAML
  are important to you**. The only exception to this rule is if there
  are any comment lines at the top of the file, before the actual code,
  these lines are preserved as is, exactly, including blank lines.
- All CWL fields are ordered in a fixed manner. The ordering of the
  fields is defined in a set of YAML files available for each
  version of the CWL standard. The YAML files are described below. 
- All strings that fit within 80 columns are expressed in flow style.
  Longer strings or strings with new lines are expressed in block style.
- All lists and maps are expressed in block fashion unless they are
  simple strings and are short enough to fit within 80 columns, in which
  case they are expressed in flow style.
- The ordering of all lists are preserved
- Maps are ordered as follows
  - Maps with no intrinsic ordering, like I/O ports, are ordered
    alphabetically
  - Steps expressed as a map are ordered by topologically sorting the
    workflow DAQ. i.e. steps that execute first are ordered first. This 
    allows step ordering to be stable over id changes.
    If steps can execute in parallel they are ordered alphabetically. 

### Inference rules

The template for a CWL version is organized as follows 

``` 
Workflow:
    - simple_field
    - complex_field
        - 
            - f11
            - f12
        -
            - f21
            - f22

CommandLineTool:
    - 

ExpressionTool:
    -
```

A simple type inference rule is used: 
