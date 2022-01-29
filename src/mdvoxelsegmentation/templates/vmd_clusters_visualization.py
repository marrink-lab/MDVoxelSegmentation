#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 12:38:38 2019

This file is used to run within VMD with python2 and numpy compiled.
After succesfully running this file you can used the 'user' selection
identiyfer in VMD followed by a cluster number to visualize that cluster.
You can save the visualization state after having the right settings in 
vmd, but you will always have to excecute this file first before loading
the visualization state (you should remove the load .gro and .xtc from
the VMD statefile).  

@author: bart
"""
from Molecule import Molecule
# needed for python2 VMD
#from atomsel import *
# needed for opening the array.npy (possibly compressed)
import numpy as np

mol = Molecule()
# The .gro file should have the same order as the tpr used for clustering.
mol.load("your_ref.gro")
mol.delFrame() # Reading in a .gro file adds coordinates, which you don't want.
mol.load("your_ref.xtc") # This should be the same .xtc used for clustering.

# The cluster file generated to visualize in vmd.
clusters = np.load('clusters.npy')
asel = atomsel("all")

for frame, cluster in enumerate(clusters):
    try:
        asel.frame = frame
    except IndexError:
        print ('Not all frames could have their clusters assigned. Stopped in \
frame ' + str(frame) + '.')
        break
    print (cluster, len(asel), len(cluster), frame)
    asel.set('user', cluster)
