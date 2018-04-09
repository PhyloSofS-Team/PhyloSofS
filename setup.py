# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name="phylosofs",
    packages=["phylosofs"],
    description="A tool to model the evolution and structural impact of "
                "alternative splicing.",
    version="0.1",
    url="https://github.com/elolaine/PhyloSofS",
    author="Adel Ait-hamlat, Lélia Polit, Diego Javier Zea, Hugues Richard "
           "and Elodie Laine",
    author_email="elodie.laine@upmc.fr",
    keywords=["splicing", "evolution", "structure"],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 2",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=[
        "numpy",      # handle transcript tables
        "networkx",   # handle phylogenetic trees
        "pydot"       # visualization of the phylogenies
    ],
    entry_points={
        'console_scripts': [
            'phylosofs=phylosofs.phylosofs:main'
        ],
    },

)
