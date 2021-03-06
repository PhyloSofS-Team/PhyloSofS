# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='phylosofs',
    description='A tool to model the evolution and structural impact of '
    'alternative splicing.',
    keywords=['splicing', 'evolution', 'structure'],
    version='0.1.0',
    url='https://github.com/PhyloSofS-Team/PhyloSofS',
    author='Adel Ait-hamlat, Diego Javier Zea, Antoine Labeeuw, Lélia Polit, '
    'Hugues Richard and Elodie Laine',
    author_email='elodie.laine@upmc.fr',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English', 'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    entry_points={
        'console_scripts': [
            'phylosofs=phylosofs.phylosofs:main',
            'setup_databases=phylosofs.setup_databases:main'
        ],
    },
    package_data={
        'phylosofs': [
            'src/plots_with_exons.jl', 'src/reconstruct_pir.jl',
            'src/reconstruct_plot.jl', 'src/setup_databases.jl',
            'src/get_pdbs.jl', 'src/Manifest.toml', 'src/Project.toml'
        ]
    },
    packages=find_packages(include=['phylosofs']),
    setup_requires=['pytest-runner'],
    install_requires=[
        'biopython',  # change B factor
        'numpy',  # handle transcript tables
        'networkx==2.3.0',  # handle phylogenetic trees
        'pydot'  # visualization of the phylogenies
    ],
    test_suite='tests',
    tests_require=[
        'pytest', 'pytest-cov', 'coveralls', 'codecov', 'pytest-pylint',
        'pylint'
    ],
    license='MIT license')
