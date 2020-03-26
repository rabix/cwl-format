# Copyright (c) 2020 Seven Bridges

import pathlib

from cwlformat.formatter import cwl_format, yaml, format_node


current_path = pathlib.Path(__file__).parent


def test_formatting_battery():

    path = current_path / "cwl"
    for raw_name in path.glob("original-*"):
        expected_name = raw_name.parent / pathlib.Path("formatted-" + "-".join(raw_name.stem.split("-")[1:]) + ".cwl")
        formatted_cwl = cwl_format(raw_name.open("r").read())
        expected_raw_cwl = expected_name.open("r").read()

        assert formatted_cwl == expected_raw_cwl


def test_node_conservation():
    path = current_path / "cwl"
    for raw_name in path.glob("original-*"):
        original_cwl = yaml.load(raw_name.open("r").read())
        formatted_cwl = format_node(original_cwl, node_path=[])

        assert formatted_cwl == original_cwl
