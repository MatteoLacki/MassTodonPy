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

from setuptools import setup, find_packages

setup(
    name          = 'MassTodonPy',
    packages      = find_packages(),
    version       = '0.3.3',
    description   = 'Estimate the products of Electron Transfer Dissociation in Mass Spectrometry for a given biological substance and the chemical reaction probabilities that lead to these products.',
    author        = 'Mateusz Krzysztof Lacki',
    author_email  = 'matteo.lacki@gmail.com',
    url           = 'https://github.com/MatteoLacki/MassTodonPy',
    download_url  = 'https://github.com/MatteoLacki/MassTodonPy/archive/0.3.3.tar.gz',
    keywords      = [
        'Mass Spectrometry',
        'ETD',
        'Electron Transfer Dissociation',
        'Fragmentation'],
    classifiers   = [
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 2.7'],
    install_requires=[
        'linearCounter',
        'numpy',
        'pandas',
        'pyteomics>=3.4.1',
        'lxml',
        'cvxopt',
        'IsoSpecPy',
        'networkx',
        'intervaltree'],
    scripts=[
        'bin/masstodon',
        'bin/masstodon_example_call'  ],
    include_package_data = True,
    package_data={
        "Data": ["Data/isotopes.txt",
                 "Data/amino_acids.txt",
                 "Data/substanceP.example",
                 "Data/ubiquitin.example"],
    }
)
