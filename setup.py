import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
	long_description = readme.read()

with open("requirements.txt", 'r') as requirements:
	install_requires = list(requirements.read().splitlines())

setup(
	name='wholecell-lens',
	version='0.0.2',
	packages=['lens'],
	author='Eran Agmon',
	author_email='eagmon@stanford.edu',
	url='https://github.com/CovertLab/Lens',
	license='MIT',
	long_description=long_description,
	long_description_content_type='text/markdown',
	install_requires=install_requires
	)
