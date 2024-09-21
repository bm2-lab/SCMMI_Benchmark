#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
requirements = [
    "numpy>=1.16.2",
    "scanpy>=1.4.6"
]

author = (
    "Shaliu Fu"
)


setup(
    author=author,
    author_email="23310100@tongji.edu.cn",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: bioinformatics",
    ],
    description="single-cell multimodal integration benchmark",
    install_requires=requirements,
    license="MIT license",
    # long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="scmmib",
    name="scmmib",
    packages=find_packages(),
    test_suite="test",
    url="https://github.com/bm2-lab/SCMMI_Benchmark",
    version="0.0.1",
    zip_safe=False,
)