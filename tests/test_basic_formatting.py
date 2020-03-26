# Copyright (c) 2020 Seven Bridges

import pathlib

from cwlformat.formatter import cwl_format


current_path = pathlib.Path(__file__).parent


def test_top_comment():
    cwl_path = pathlib.Path(current_path / "cwl" / "original-commented.cwl")
    raw_cwl = open(cwl_path, "r").read()
    formatted_cwl = cwl_format(raw_cwl)

    first_lines_raw = raw_cwl.splitlines(keepends=True)[:3]
    first_lines_formatted = formatted_cwl.splitlines(keepends=True)[:3]

    assert first_lines_raw == first_lines_formatted


def test_node_conservation():
    pass
