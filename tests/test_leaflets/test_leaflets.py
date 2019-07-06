#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 14:35:52 2019

@author: bart
"""
import itertools
import numpy as np
import MDAnalysis as mda
from mdvoxelclustering import leaflets


def test_contour_clustering_bilayer():
    """
    Tests the inner contour clustering of an atomgroup for a single bilayer, 
    without the presence of an exclusion mask.
    """
    ref_cluster1 = set(itertools.product(list(range(10)), list(range(10)), [6]))
    ref_cluster2 = set(itertools.product(list(range(10)), list(range(10)), [3]))
    
    universe = mda.Universe('tests/test_leaflets/test_bilayer.tpr', 
                            'tests/test_leaflets/test_bilayer.gro')
    lipids = universe.select_atoms('resname DOPC')
    output = leaflets.contour_clustering(
        lipids, exclusion_mask=False, resolution=1, span=0,
        inv=True)
    clusters = output[0]
    
    assert (set(clusters.get(1)) == ref_cluster1 and
            set(clusters.get(2)) == ref_cluster2 and
            len(clusters.get(1)) == len(set(clusters.get(1))) and
            len(clusters.get(2)) == len(set(clusters.get(2)))
            )


def test_volume_clustering_bilayer():
    """
    Tests the volume clustering of an atomgroup for a single bilayer, 
    without the presence of an exclusion mask.
    """
    ref_cluster1 = set(itertools.product(list(range(10)), list(range(10)), 
                                         list(range(3, 7))))
    
    universe = mda.Universe('tests/test_leaflets/test_bilayer.tpr', 
                            'tests/test_leaflets/test_bilayer.gro')
    lipids = universe.select_atoms('resname DOPC')
    output = leaflets.volume_clustering(
            lipids, headgroups_selection = False,
            exclusion_mask = False, resolution = 1
            )
    clusters = output[0]
    
    assert (set(clusters.get(1)) == ref_cluster1 and
            len(clusters.get(1)) == len(set(clusters.get(1)))
            )
    
    
def test_leaflet_clustering():
    """
    Testing the leaflet clustering on a single bilayer of only DOPC.
    """
    universe = mda.Universe('tests/test_leaflets/test_bilayer.tpr', 
                            'tests/test_leaflets/test_bilayer.gro')
    #universe = mda.Universe('test_bilayer.tpr', 
    #                        'test_bilayer.gro')
    lipids = universe.select_atoms('resname DOPC')
    tails_selection = universe.select_atoms('name C3A C4A C5A D3A D4A D5A '
                                            'C3B C4B C5B D3B D4B D5B')
    headgroups_selection = universe.select_atoms('name PO4 NH3 GL1 GL2')
    
    ref_leaflet1 = universe.residues[0:169].atoms.ix
    ref_leaflet2 = universe.residues[169:338].atoms.ix
    ref_leaflets = np.zeros(len(lipids.universe.atoms))
    ref_leaflets[ref_leaflet1] = 1
    ref_leaflets[ref_leaflet2] = 2
    
    clusters = leaflets.leaflet_clustering(
        tails_selection, headgroups_selection,
        exclusions_selection = False, 
        resolution = 1, bits = 'uint32', verbose = False,
        )
    
    assert np.all(clusters == ref_leaflets)
    
    
