import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

with open("requirements.txt", 'r') as requirements:
    install_requires = list(requirements.read().splitlines())

setup(
    name='vivarium',
    version='0.0.28',
    packages=[
        'vivarium',
        'vivarium.actor',
        'vivarium.analysis',
        'vivarium.environment',
        'vivarium.processes',
        'vivarium.data',
        'vivarium.data.flat',
        'vivarium.data.flat.media',
        'vivarium.data.json_files',
        'vivarium.surrogates',
        'vivarium.utils',
        'vivarium.utils._netflow'],
    author='Eran Agmon, Ryan Spangler',
    author_email='eagmon@stanford.edu, spanglry@stanford.edu',
    url='https://github.com/CovertLab/Lens',
    license='MIT',
    entry_points={
        'console_scripts': [
            'vivarium.actor.boot=vivarium.actor.boot:run',
            'vivarium.actor.control=vivarium.actor.control:run',
            'vivarium.environment.boot=vivarium.environment.boot:run',
            'vivarium.environment.control=vivarium.environment.control:run',
            'vivarium.composites=vivarium.composites:run']},
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={
        'vivarium.data.flat': ['*.tsv'],
        'vivarium.data.flat.media': ['*.tsv'],
        'vivarium.data.json_files': ['*.json']},
    include_package_data=True,
    install_requires=install_requires)
