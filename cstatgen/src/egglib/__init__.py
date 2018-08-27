"""
EggLib
======

The Python version (egglib-py) provides access to classes of the
EggLib C++ library (through a SWIG wrapping) and an extended
collection of Python tools through a simplified interface.

"""


__license__ = """
    Copyright 2008-2013 Stephane De Mita, Mathieu Siol

    This file is part of EggLib.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
"""

import utils
from data import (SequenceItem, Container, Align, GenBank, GenBankFeature,
                  GenBankFeatureLocation, SSR, TIGR, Tree, TreeNode)



version = '2.1.7'
""" egglib version number. """

