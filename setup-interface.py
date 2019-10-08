# This file provides a stripped down version of lens that declares only
# the interface defined in the lens.actor package, so that clients
# may implement the lens interface without having to be dependent on
# the rest of the lens package. 

import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

with open("requirements-interface.txt", 'r') as requirements:
    install_requires = list(requirements.read().splitlines())

setup(
    name='lens-interface',
    version='0.0.27',
    packages=[
        'lens',
        'lens.actor'],
    author='Eran Agmon, Ryan Spangler',
    author_email='eagmon@stanford.edu, spanglry@stanford.edu',
    url='https://github.com/CovertLab/Lens',
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires)
