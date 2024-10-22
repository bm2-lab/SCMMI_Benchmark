#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
requirements = [
    "numpy>=1.16.2",
    "scanpy>=1.4.6",
    "scib",
    "scanpy",
    "scglue"
]


extras_requirements = {
    "notebooks": [
        "louvain>=0.6.1",
        "python-igraph>=0.7.1.post6",
        "colour>=0.1.5",
        "umap-learn>=0.3.8",
        "seaborn>=0.9.0",
        "leidenalg>=0.7.0",
    ],
    "docs": [
        "sphinx>=1.7.1",
        "nbsphinx",
        "sphinx_autodoc_typehints",
        "sphinx-rtd-theme",
        "myst_parser",
        "sphinx-automodapi",
    ],
}

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
    description="Single-cell multimodal integration benchmark",
    install_requires=requirements+extras_requirements['docs']+extras_requirements['notebooks'],
    license="MIT license",
    # long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="scmmib",
    name="scmmib",
    packages=find_packages(),
    test_suite="test",
    extras_require=extras_requirements,
    url="https://github.com/bm2-lab/SCMMI_Benchmark",
    version="0.0.1",
    zip_safe=False,
)