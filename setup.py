from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
        long_description = fh.read()

setup(name='selfcalframework',
version='0.1.3',
url='https://github.com/miguelcarcamov/objectoriented_selfcal',
description='A Python object oriented framework to do self-calibration',
author='Miguel Carcamo',
author_email='miguel.carcamo@manchester.ac.uk',
long_description=long_description,
long_description_content_type="text/markdown",
license='GNU',
packages=find_packages(),
classifiers=[
"Programming Language :: Python :: 3",
"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
"Operating System :: OS Independent"],
)
