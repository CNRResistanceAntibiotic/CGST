#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the setup file of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = ""
exec(open('cgst/version.py').read())


setup(name='CGST',
      version=__version__,
      description='CGST: Core-Genome Sequence Typing',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/CNRResistanceAntibiotic/CGST',
      author='Aurélien Birer',
      author_email='abirer36@gmail.com',
      license='GPLv3',
      packages=["cgst"],
      include_package_data=True,
      install_requires=['Bio'],
      entry_points={"console_scripts": ['cgst = cgst.CGST:run']},
      zip_safe=False,
      python_requires='>=3.6')
