import os
import setuptools
import sys
import platform

from pybind11.setup_helpers import Pybind11Extension, build_ext

extra_compile_args = []
if sys.platform == 'win32':
  extra_compile_args += [
    '/std:c++20', '/O2'
  ]
else:
  extra_compile_args += [
    '-std=c++20', '-O2',
  ]

setuptools.setup(
  setup_requires=['pbr','pybind11','numpy'],
  cmdclass={"build_ext": build_ext},
  ext_modules=[
    Pybind11Extension(
        "fastosteoid",
        ["src/fastosteoid.cpp"],
        extra_compile_args=extra_compile_args,
        language="c++",
    ),
  ]
)
