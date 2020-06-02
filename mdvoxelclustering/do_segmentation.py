#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:17:19 2020

@author: bart
Running this file is all that is needed to perform a basic leaflet
segmentation for Martini lipids.
"""
from . import argparser
from . import leaflets

def main():
    # Reading in the terminal commands/input files
    args = argparser.main()

if __name__=='__main__':
    main()
