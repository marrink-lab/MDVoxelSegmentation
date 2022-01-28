#!/usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

setuptools.setup(
    packages=['mdvoxelsegmentation'],
    entry_points = {
        'console_scripts' : ['mdvseg=mdvoxelsegmentation.do_segmentation:main'],
        },
)
