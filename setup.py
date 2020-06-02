#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt', 'r') as requirements_file:
	requirements = [x for x in requirements_file.readlines()]
#test_requirements = [str(tr.req) for tr in parsed_test_requirements]


setup(
    name='mdvoxelclustering',
    version='0.9.2',
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
    license="Apache 2",
    zip_safe=False,
    keywords='mdvoxelclustering',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache 2',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    test_suite='tests',
#    tests_require=test_requirements
)
