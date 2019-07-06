#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pip

from pip._internal.req import parse_requirements

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

parsed_requirements = parse_requirements(
    'requirements/prod.txt',
    session=pip._internal.download.PipSession()
)

parsed_test_requirements = parse_requirements(
    'requirements/test.txt',
    session=pip._internal.download.PipSession()
)


requirements = [str(ir.req) for ir in parsed_requirements]
test_requirements = [str(tr.req) for tr in parsed_test_requirements]


setup(
    name='mdvoxelclustering',
    version='0.1.0',
    description="Using neighbour clustering in voxelspace for fast and consistant spatial and temporal clustering.",
    long_description=readme + '\n\n' + history,
    author="Bart M. H. Bruininks",
    author_email='b.m.h.bruininks@rug.nl',
    url='https://github.com/BartBruininks/mdvoxelclustering',
    packages=[
        'mdvoxelclustering',
    ],
    package_dir={'mdvoxelclustering':
                 'mdvoxelclustering'},
    include_package_data=True,
    install_requires=requirements,
    license="ISCL",
    zip_safe=False,
    keywords='mdvoxelclustering',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
