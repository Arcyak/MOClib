from setuptools import setup, find_packages
from __init__ import __version__

setup(
    author='Yannick HÉNIN',
    name='MOClib',
    github='https://github.com/Arcyak',
    description='Python custom MOC library that creates MOC objects from HEALPix',
    version=__version__,
    keywords=['HEALPix', 'MOC', 'Astropy'],
    packages=find_packages(),
    python_requires='>=3.7'
)
