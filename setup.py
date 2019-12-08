#!/usr/bin/env python3
from distutils.core import setup, Extension

# chemm
setup(name='compchem',
      version='0.3',
      packages=['compchem', 'compchem/view'],
      scripts=['compchem/pdb.py',
               'compchem/view/trackball.py'],
      author='drccr',
      author_email='drccr@pm.me',
      url='https://github.com/drccr/compchem-tools')
