# This file provides a stripped down version of vivarium that declares only
# the interface defined in the vivarium.actor package, so that clients
# may implement the vivarium interface without having to be dependent on
# the rest of the vivarium package.

import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

with open("requirements-interface.txt", 'r') as requirements:
    install_requires = list(requirements.read().splitlines())

setup(
    name='vivarium-interface',
    version='0.0.28',
    packages=[
        'vivarium',
        'vivarium.actor'],
    author='Eran Agmon, Ryan Spangler',
    author_email='eagmon@stanford.edu, spanglry@stanford.edu',
    url='https://github.com/CovertLab/Lens',
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires)
