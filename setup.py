#!/usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

setuptools.setup(
    packages=['mdvoxelsegmentation'],
    include_package_data=True,
    package_data={'':['templates/*']},
    entry_points={
        'console_scripts' : ['mdvseg=mdvoxelsegmentation.do_segmentation:main'],
        },
)
