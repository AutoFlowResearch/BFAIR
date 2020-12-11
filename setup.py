# -*- coding: utf-8 -*-
import re
from setuptools import setup
import setuptools

with open("README.md", "r") as doc:
    long_description = doc.read()


version = re.search(
    r'^__version__\s*=\s*"(.*)"',
    open("AutoFlow_OmicsDataHandling/AutoFlow_OmicsDataHandling.py").read(),
    re.M,
).group(1)


setup(
    name="AutoFlow_OmicsDataHandling",
    # version=version,
    author="AutoFlow",
    author_email="TBD",
    description=(
        "General Repository for Omics Data Handling tools"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Biosustain/AutoFlow-OmicsDataHandling",
    packages=setuptools.find_packages(),
    # entry_points={
    #     "console_scripts": [
    #         "AutoFlow_OmicsDataHandling =
    # AutoFlow_OmicsDataHandling.__main__:main"
    #     ]
    # },
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "flake8",
        "pytfa",
        "pyopenms",
    ],
    # dependency_links=
    include_package_data=True,
)
