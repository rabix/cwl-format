# Copyright (c) 2020 Seven Bridges

from typing import Union
import sys
import importlib.resources as pkg_resources

import ruamel.yaml
from ruamel.yaml.compat import StringIO

from cwlformat.version import __version__

yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=2, offset=0)
Literal = ruamel.yaml.scalarstring.LiteralScalarString

key_order_dict = yaml.load(pkg_resources.read_text("cwlformat", "keyorder.yml"))


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


def format_node(cwl: Union[dict, list, str], node_path=None):
    if isinstance(cwl, str):
        if len(cwl) > 80:
            return Literal(cwl)
        else:
            return cwl

    elif isinstance(cwl, dict):
        return {k: format_node(v, node_path + [k]) for k, v in reorder_node(cwl, node_path)}

    elif isinstance(cwl, list):
        return [format_node(v, node_path) for v in cwl]

    else:
        return cwl


def reorder_node(cwl: dict, node_path: list) -> dict:
    known_key_order = key_order_dict.get(
        infer_type(cwl, node_path), key_order_dict["generic-ordering"])
    extra_keys = sorted(set(cwl.keys()) - set(known_key_order))

    for k in known_key_order + extra_keys:
        if k in cwl:
            yield k, cwl[k]


def infer_type(cwl: dict, node_path: list):
    if "class" in cwl:
        return cwl["class"]
    else:
        return "generic-ordering"


def cwl_format(raw_cwl: str) -> str:
    as_dict = yaml.load(raw_cwl)
    as_dict = format_node(as_dict, node_path=[])
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
    sys.stdout.write(cwl_format(open(args.cwlfile, "r").read()))


if __name__ == "__main__":
    main()
