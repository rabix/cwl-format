# Copyright (c) 2020 Seven Bridges

import pathlib

from cwlformat.formatter import cwl_format


current_path = pathlib.Path(__file__).parent


def test_formatting_battery():

    path = current_path / "cwl"
    for raw_name in path.glob("original-*"):
        expected_name = raw_name.parent / pathlib.Path("formatted-" + "-".join(raw_name.stem.split("-")[1:]) + ".cwl")
        formatted_cwl = cwl_format(raw_name.open("r").read())
        expected_raw_cwl = expected_name.open("r").read()

        assert formatted_cwl == expected_raw_cwl
