# -*- coding: utf-8 -*-

import pathlib
from setuptools import setup, find_packages

current_path = pathlib.Path(__file__).parent
ver_path = pathlib.Path(current_path, "cwlformat", "version.py")
_ver = {}
exec(ver_path.open("r").read(), _ver)
version = _ver["__version__"]

readme = pathlib.Path(current_path, "Readme.md").read_text()

setup(
    name='cwlformat',
    python_requires='>=3.7.0',
    version=version,
    description='A prettifier for CWL code',
    long_description=readme,
    long_description_content_type="text/markdown",
    author='Kaushik Ghose',
    author_email='kaushik.ghose@sbgenomics.com',
    url='https://github.com/kaushik-work/cwlformat',
    packages=find_packages(exclude=('tests', 'docs')),
    entry_points={
        'console_scripts': [
            'cwl-format=cwlformat.formatter:main'
        ],
    },
    include_package_data=True,
    install_requires=[
        "ruamel.yaml >= 0.15.77",
    ]
)
