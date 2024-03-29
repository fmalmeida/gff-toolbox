#!/usr/bin/env python3
"""
This script allows a user to run gff-toolbox without a pip/setuptools installation. Assuming all
dependencies are installed, they can simply clone the repo and run this script.

Copyright 2020 Felipe Almeida (almeidafmarques@gmail.com)
https://github.com/fmalmeida/gff-toolbox

This file is part of gff-toolbox. gff-toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. gff-toolbox is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with gff-toolbox.
If not, see <http://www.gnu.org/licenses/>.
"""

from gfftoolbox.__main__ import main


if __name__ == '__main__':
    main()
