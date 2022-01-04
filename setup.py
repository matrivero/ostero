from setuptools import find_packages

from numpy.distutils.core import setup, Extension

ext1 = Extension(name='external.modf90',
                 sources=['src/external.f90'],
                 f2py_options=['--quiet'],
                )


setup(name="ostero",
      version="0.0.1",
      package_dir={"": "src"},
      packages=find_packages(where="src"),
      ext_modules=[ext1])
