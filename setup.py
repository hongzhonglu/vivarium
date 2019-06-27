import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
	long_description = readme.read()

current_dir = os.getcwd()
lens_dir = os.path.join(current_dir, 'lens')

setup(
	name='lens',
	version='0.0.1',
	packages=['lens'],
	author='Eran Agmon',
	author_email='eagmon@stanford.edu',
	url='https://github.com/CovertLab/Lens',
	license='MIT',
	long_description=long_description,
	long_description_content_type='text/markdown')
