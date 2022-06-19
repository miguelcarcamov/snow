import distutils.text_file
# read the contents of your README file
from pathlib import Path
from typing import List

from setuptools import find_packages, setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


def _parse_requirements(filename: str) -> List[str]:
    """Return requirements from requirements file."""
    # Ref: https://stackoverflow.com/a/42033122/
    return distutils.text_file.TextFile(
        filename=str(Path(__file__).with_name(filename))).readlines()


setup(
    name='selfcalframework',
    version='0.1.3',
    url='https://github.com/miguelcarcamov/objectoriented_selfcal',
    description='A Python object oriented framework to do self-calibration',
    author='Miguel Carcamo',
    author_email='miguel.carcamo@manchester.ac.uk',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=_parse_requirements("requirements.txt"),
    license='GNU',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
)
