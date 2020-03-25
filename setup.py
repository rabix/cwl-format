# -*- coding: utf-8 -*-

import pathlib
from setuptools import setup, find_packages


with open('Readme.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license_text = f.read()

current_path = pathlib.Path(__file__).parent
ver_path = pathlib.Path(current_path, "cwlformat", "version.py")
_ver = {}
exec(ver_path.open("r").read(), _ver)
version = _ver["__version__"]


setup(
    name='CWLformat',
    python_requires='>=3.7.0',
    version=version,
    description='A prettifier for CWL code',
    long_description=readme,
    author='Kaushik Ghose',
    author_email='kaushik.ghose@sbgenomics.com',
    url='https://github.com/kaushik-work/cwlformat',
    license=license_text,
    packages=find_packages(exclude=('tests', 'docs')),
    entry_points={
        'console_scripts': [
            'cwl-format=cwlformat.formatter:main'
        ],
    },
    # https://docs.python.org/3.7/distutils/setupscript.html#installing-additional-files
    data_files=[("keyorder", ["cwlformat/keyorder.yml"])],
    install_requires=[
        "ruamel.yaml >= 0.15.77",
    ]
)
