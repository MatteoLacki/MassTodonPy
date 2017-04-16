# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

# from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name          = 'MassTodonPy',
    # packages      = ['MassTodonPy', 'MassTodonPy.data', 'MassTodonPy.Deconvolutor'], # this must be the same as the name above
    packages      = find_packages(),
    version       = '0.21',
    description   = 'Estimate the products of Electron Transfer Dissociation in Mass Spectrometry for a given biological substance and the chemical reaction probabilities that lead to these products.',
    author        = 'Mateusz Krzysztof Lacki',
    author_email  = 'matteo.lacki@gmail.com',
    url           = 'https://github.com/MatteoLacki/MassTodonPy',
    download_url  = 'https://github.com/MatteoLacki/MassTodonPy/archive/0.21.tar.gz',
    keywords      = ['Mass Spectrometry', 'ETD', 'Electron Transfer Dissociation', 'Fragmentation'],
    classifiers   = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 2.7'
    ],
    install_requires=[
        'linearCounter', 'numpy', 'pandas', 'pyteomics>=3.4.1', 'cvxopt', 'IsoSpecPy', 'scipy', 'matplotlib', 'python-igraph', 'networkx', 'intervaltree'
    ],
    include_package_data = True,
    package_data={
        "MassTodonPy": [
            "data/isotopes.txt"
        ],
    }
)
