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
    version='0.0.12',
    packages=[
        'lens',
        'lens.actor',
        'lens.analysis',
        'lens.environment',
        'lens.environment.condition',
        'lens.environment.condition.look_up_tables',
        'lens.environment.condition.look_up_tables.transport_concentrations',
        'lens.environment.condition.look_up_tables.transport_fluxes',
        'lens.environment.condition.media',
        'lens.processes',
        'lens.reconstruction',
        'lens.reconstruction.CovertPalsson2002',
        'lens.reconstruction.flat',
        'lens.reconstruction.kinetic_rate_laws',
        'lens.reconstruction.kinetic_rate_laws.parameters',
        'lens.surrogates',
        'lens.utils',
        'lens.utils._netflow',
        'lens.utils.io'],
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
        'lens.environment.condition': ['*.tsv'],
        'lens.environment.condition.media': ['*.tsv'],
        'lens.environment.condition.look_up_tables.transport_concentrations': ['*.tsv'],
        'lens.environment.condition.look_up_tables.transport_fluxes': ['*.tsv'],
        'lens.reconstruction.CovertPalsson2002': ['*.tsv'],
        'lens.reconstruction.flat': ['*.tsv'],
        'lens.reconstruction.kinetic_rate_laws.parameters': ['*.json']},
    include_package_data=True,
    install_requires=install_requires)
