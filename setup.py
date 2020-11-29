# Setup file for the INCA_MFA_tools
# -*- coding: utf-8 -*-
import re
from setuptools import setup
import setuptools

with open("README.md", "r") as doc:
    long_description = doc.read()

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('INCA_script_generator/INCA_script_generator.py').read(),
    re.M
    ).group(1)

setup(
    name = "INCA_script_generator",
    version = version,
    author = "Matthias Mattanovich",
    author_email = "matmat@biosustain.com",
    description = ("Generates an MFA script for INCA, a MATLAB software, that also calls MATLAB"),
    long_description = long_description, 
    long_description_content_type="text/markdown",
    url = "https://github.com/mmattano/INCA_script_generator",
    packages=setuptools.find_packages(),
    entry_points = {
        "console_scripts": ["INCA_script_generator = INCA_script_generator.__main__:main"]
        },
    license = "MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
    install_requires=['pandas', 'numpy'],
    include_package_data = True,

)