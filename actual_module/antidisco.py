#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:35:54 2019

@author: bart
"""
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import MDAnalysis as mda

def plot_clusters(universe, clusters, skip, reduce_points):
    try:
        os.mkdir('figs')
    except FileExistsError:
        pass
    print('Some basic plotting...')
    for frame_idx, _ in enumerate(universe.trajectory[::skip]):
        frame = universe.trajectory.frame
        fig = plt.figure(figsize = [10, 10])
        fig.suptitle('Frame {} at {} ns, with {} clusters.'.format(
                frame, 
                frame*universe.trajectory.dt/1000, 
                len(set(clusters[frame_idx]))-1),
                    )
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        ax.set_xlim3d(0, universe.dimensions[0])
        ax.set_ylim3d(0, universe.dimensions[1])
        ax.set_zlim3d(0, universe.dimensions[2])
        #print('Current frame: {}'.format(frame_idx))
        for cluster in set(clusters[frame_idx]):
            if cluster == 0:
                continue
            #print(len(clusters[frame_idx][clusters[frame_idx] == cluster]))
            ax.scatter(universe.atoms.positions[clusters[frame_idx] == cluster][:,0][::reduce_points], universe.atoms.positions[clusters[frame_idx] == cluster][:,1][::reduce_points], universe.atoms.positions[clusters[frame_idx] == cluster][:,2][::reduce_points], alpha = 0.1)
        fig.savefig('figs/leaflets_frame-{:09d}.png'.format(universe.trajectory.frame), dpi = 300)
        plt.close()
        
if __name__=='__main__':
    print('Loading universe...')
    data = mda.Universe('/coarse/bart/projects/clustering/test_files/4_transfection/md.tpr', '/coarse/bart/projects/clustering/test_files/4_transfection/centered_trajout.xtc')
    print('Loading clusters...')
    clusters = np.load('/home/bart/projects/MDVoxelClustering/clusters_ordered.npy')
    plot_clusters(data, clusters, 1, 10)
        
