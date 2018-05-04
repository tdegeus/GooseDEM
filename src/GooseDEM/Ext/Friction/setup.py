desc = '''
GooseDEM is a C++ module, wrapped in Python, that provides several simple methods to run a (small)
DEM simulation.
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import cppmat

header = open('Friction.h','r').read()
world  = re.split('(.*)(\#define GOOSEDEM_EXT_FRICTION_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split('(.*)(\#define GOOSEDEM_EXT_FRICTION_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split('(.*)(\#define GOOSEDEM_EXT_FRICTION_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GooseDEM_Ext_Friction',
    ['python.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True ),
      cppmat  .find_eigen()
    ],
    language='c++'
  ),
]

setup(
  name             = 'GooseDEM_Ext_Friction',
  description      = 'DEM simulations',
  long_description = desc,
  version          = __version__,
  license          = 'MIT',
  author           = 'Tom de Geus, Tianxiao Ma, Wencheng Ji',
  author_email     = 'tom@geus.me',
  url              = 'https://github.com/tdegeus/GooseDEM',
  ext_modules      = ext_modules,
  install_requires = ['pybind11>=2.2.0','cppmat>=0.4.1'],
  cmdclass         = {'build_ext': cppmat.BuildExt},
  zip_safe         = False,
)
