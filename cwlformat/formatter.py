# Copyright (c) 2020 Seven Bridges

from typing import Union

import ruamel.yaml
from ruamel.yaml.compat import StringIO

from cwlformat.version import __version__

yaml = ruamel.yaml.YAML()
Literal = ruamel.yaml.scalarstring.LiteralScalarString


def leading_comment_lines(raw_cwl: str):

    if len(raw_cwl) > 0 and raw_cwl.lstrip()[0] == "{":
        return ""

    top_comment = []
    for _line in raw_cwl.splitlines(keepends=True):
        line = _line.strip()
        if line == "" or line[0] == "#":
            top_comment += [_line]
        else:
            break

    return "".join(top_comment)


def format_node(cwl: Union[dict, list, str]):
    if isinstance(cwl, str):
        if len(cwl) > 80:
            return Literal(cwl)
        else:
            return cwl

    elif isinstance(cwl, dict):
        return {k: format_node(v) for k, v in reorder_node(cwl)}

    elif isinstance(cwl, list):
        return [format_node(v) for v in cwl]

    else:
        return cwl


def reorder_node(cwl: dict) -> dict:
    for k in sorted(cwl.keys()):
        yield k, cwl[k]


def cwl_format(raw_cwl: str) -> str:
    as_dict = yaml.load(raw_cwl)
    as_dict = format_node(as_dict)
    stream = StringIO()
    yaml.dump(as_dict, stream)
    return leading_comment_lines(raw_cwl) + stream.getvalue()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"Rabix/cwl-format v{__version__}\n"
                    "A very opinionated code formatter for CWL")
    parser.add_argument("cwlfile")
    args = parser.parse_args()
    print(cwl_format(open(args.cwlfile, "r").read()))


if __name__ == "__main__":
    main()
