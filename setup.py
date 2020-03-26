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
    name='cwlformat',
    python_requires='>=3.7.0',
    version=version,
    description='A prettifier for CWL code',
    long_description=readme,
    long_description_content_type="text/markdown",
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
    include_package_data=True,
    install_requires=[
        "ruamel.yaml >= 0.15.77",
    ]
)
