# Copyright (c) 2020 Seven Bridges

import pathlib

import ruamel.yaml

from cwlformat.explode import explode, CWLProcess

yaml = ruamel.yaml.YAML()

current_path = pathlib.Path(__file__).parent


def test_explode():
    path = current_path / "cwl"
    src_fp = path / "formatted-atac-seq-pipeline.cwl"
    fp_out = src_fp.parent / "expected-exploded-atac-seq.cwl"

    as_dict = yaml.load(src_fp.read_text())
    for exploded in explode(CWLProcess(as_dict, fp_out)):
        assert exploded.file_path.exists()
