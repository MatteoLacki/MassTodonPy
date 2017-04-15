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

# I have no idea what I am doing...
setup(
    name          = 'MassTodonPy',
    # packages      = ['MassTodonPy', 'MassTodonPy.data', 'MassTodonPy.Deconvolutor'], # this must be the same as the name above
    packages      = ['misc'].extend(find_packages()),
    version       = '0.21',
    description   = 'A module that investigates the products of Electron Transfer Dissociation in Mass Spectrometry for a given biological substance.',
    author        = 'Mateusz Lacki',
    author_email  = 'matteo.lacki@gmail.com',
    url           = 'https://github.com/MatteoLacki/MassTodonPy',
    download_url  = 'https://github.com/MatteoLacki/MassTodonPy/archive/0.21.tar.gz',
    keywords      = ['Mass Spectrometry', 'ETD', 'Electron Transfer Dissociation', 'Fragmentation'],
    classifiers   = [],
)
