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
    version='0.0.23',
    packages=[
        'lens',
        'lens.actor',
        'lens.analysis',
        'lens.environment',
        'lens.processes',
        'lens.data',
        'lens.data.flat',
        'lens.data.flat.media',
        'lens.data.json_files',
        'lens.surrogates',
        'lens.utils',
        'lens.utils._netflow'],
    author='Eran Agmon, Ryan Spangler',
    author_email='eagmon@stanford.edu, spanglry@stanford.edu',
    url='https://github.com/CovertLab/Lens',
    license='MIT',
    entry_points={
        'console_scripts': [
            'lens.actor.boot=lens.actor.boot:run',
            'lens.actor.control=lens.actor.control:run',
            'lens.environment.boot=lens.environment.boot:run',
            'lens.environment.control=lens.environment.control:run',
            'lens.composites=lens.composites:run']},
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={
        'lens.data.flat': ['*.tsv'],
        'lens.data.flat.media': ['*.tsv'],
        'lens.data.json_files': ['*.json']},
    include_package_data=True,
    install_requires=install_requires)
