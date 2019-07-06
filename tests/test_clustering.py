#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_mdvoxelclustering
----------------------------------

Tests for `mdvoxelclustering` module.
"""
#from mdvoxelclustering import Mdvoxelclustering
import numpy as np
import MDAnalysis as mda
import mdvoxelclustering as mdv
from  mdvoxelclustering import clustering as clus

def test_forward_and_backward():
    """
    Testing forward and backward mapping.
    """
    universe = mda.Universe('test.gro')
    lipids = universe.select_atoms('resname DLPC')
    forward_atoms = lipids.atoms 
    density_mask, voxel2atom = clus.gen_explicit_matrix(lipids, resolution = 1, 
                                                        PBC = 'cubic')
    voxel_list = np.asarray(np.where(density_mask == True)).T
    backward_atoms = clus.convert_voxels2atomgroup(voxel_list, voxel2atom, lipids)
    
    assert np.all(forward_atoms.ix == backward_atoms.ix)
    
    
def test_clustering_and_backward_no_exclusions():
    """
    Testing clustering and backward mapping (no eclusions).
    """
    universe = mda.Universe('test.gro')
    lipids = universe.select_atoms('resname DLPC')
    
    density_mask, voxel2atom = clus.gen_explicit_matrix(lipids, resolution = 1, 
                                                        PBC = 'cubic')
    clusters = clus.set_clustering(density_mask)
    cluster_atomgroups = clus.convert_clusters2atomgroups(clusters, voxel2atom, 
                                                          lipids)
    cluster0_ref = set([17,18,19])
    cluster1_ref = set([5,6,7,8,9,10,11,12])
    cluster0 = set(cluster_atomgroups[0].ix)
    cluster1 = set(cluster_atomgroups[1].ix)
    
    assert (np.all(cluster0_ref == cluster0) and 
            np.all(cluster1_ref == cluster1))


def test_clustering_and_backward_exclusions():
    """
    Testing clustering and backward mapping (eclusions).
    """
    universe = mda.Universe('test.gro')
    lipids = universe.select_atoms('resname DLPC')
    exclusions = universe.select_atoms('name PO4')
    density_mask, voxel2atom = clus.gen_explicit_matrix(lipids, resolution = 1, 
                                                        PBC = 'cubic')
    density_mask_exclusions, _ = clus.gen_explicit_matrix(exclusions, 
                                                          resolution = 1, 
                                                          PBC = 'cubic')
    clusters = clus.set_clustering(density_mask, density_mask_exclusions)
    cluster_atomgroups = clus.convert_clusters2atomgroups(clusters, voxel2atom, 
                                                          lipids)
    
    cluster0_ref = set([17,18,19])
    cluster1_ref = set([12])
    cluster2_ref = set([5,6,7,8,9,10])
    cluster0 = set(cluster_atomgroups[0].ix)
    cluster1 = set(cluster_atomgroups[1].ix)
    cluster2 = set(cluster_atomgroups[2].ix)
    
    assert (np.all(cluster0_ref == cluster0) and 
            np.all(cluster1_ref == cluster1) and 
            np.all(cluster2_ref == cluster2))


def test_matrix_blurring_span0():
    """
    Testing matrix blurring with a line blur of span 1. This is a
    special case.
    """
    ref_matrix = np.zeros(27)
    ref_matrix[4] = True
    ref_matrix[10] = True
    ref_matrix[12:15] = True
    ref_matrix[16] = True
    ref_matrix[22] = True
    ref_matrix = ref_matrix.reshape([3,3,3])
    
    test_matrix = np.zeros(27)
    test_matrix[13] = True
    test_matrix = test_matrix.reshape([3,3,3])
    blurred_matrix = clus.blur_matrix(test_matrix, span = 0, PBC = 'cubic')
    
    assert np.all(blurred_matrix == ref_matrix)


def test_matrix_blurring_span_negative_one():
    """
    Testing matrix blurring with a positive line blur of span 1. This is a
    special case.
    """
    ref_matrix = np.zeros(27)
    ref_matrix[13] = True
    ref_matrix[14] = True
    ref_matrix[16] = True
    ref_matrix[22] = True
    ref_matrix = ref_matrix.reshape([3,3,3])
    
    test_matrix = np.zeros(27)
    test_matrix[13] = True
    test_matrix = test_matrix.reshape([3,3,3])
    blurred_matrix = clus.blur_matrix(test_matrix, span = -1, PBC = 'cubic')
    
    assert np.all(ref_matrix == blurred_matrix)


def test_matrix_blurring_spanN():
    """
    Testing matrix block blurring of span 1. This the general case.
    """
    ref_matrix = np.zeros(125)
    ref_matrix = ref_matrix.reshape([5,5,5])
    ref_matrix[1:4, 1:4, 1:4] = True
    
    test_matrix = np.zeros(125)
    test_matrix[62] = True
    test_matrix = test_matrix.reshape([5,5,5])
    blurred_matrix = clus.blur_matrix(test_matrix, span = 1, PBC = 'cubic')
    
    assert np.all(blurred_matrix == ref_matrix)
    

def test_matrix_blurring_cubicPBC():
    """
    Testing matrix blurring over periodic blurring with span 1 line blur.
    """
    ref_matrix = np.zeros(27)
    ref_matrix = ref_matrix.reshape([3,3,3])
    ref_matrix[0,0,0] = True
    ref_matrix[1,0,0] = True
    ref_matrix[2,0,0] = True
    ref_matrix[0,1,0] = True
    ref_matrix[0,2,0] = True
    ref_matrix[0,0,1] = True
    ref_matrix[0,0,2] = True
    
    test_matrix = np.zeros(27)
    test_matrix[0] = True
    test_matrix = test_matrix.reshape([3,3,3])
    blurred_matrix = clus.blur_matrix(test_matrix, span = 0, PBC = 'cubic')
    
    assert np.all(blurred_matrix == ref_matrix)
