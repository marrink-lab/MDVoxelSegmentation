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


setup(
    name='mdvoxelsegmentation',
    version='1.1.0',
    description="Using neighbors in voxelspace for fast and consistant spatial and temporal segmentation.",
    long_description=readme + '\n\n' + history,
    author="Bart M. H. Bruininks",
    author_email='b.m.h.bruininks@rug.nl',
    url='https://github.com/BartBruininks/MDVoxelSegmentation',
    packages=[
        'mdvoxelsegmentation',
    ],
    package_dir={'mdvoxelsegmentation':
                 'mdvoxelsegmentation'},
    include_package_data=True,
    install_requires=requirements,
    license="Apache 2",
    zip_safe=False,
    keywords=['md', 'voxel', 'segmentation', 'analysis'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    test_suite='tests',
)
