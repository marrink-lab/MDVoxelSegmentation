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
import sys

def plot_clusters(universe, clusters, skip, reduce_points, colors, start_frame = 0, stop_frame = False, min_cluster_size = 150):
    try:
        os.mkdir('figs_sorted')
    except FileExistsError:
        pass
    print('Some basic plotting...')
    if stop_frame == False:
        stop_frame = len(universe.trajectory)
    for frame_idx, _ in enumerate(universe.trajectory[start_frame:stop_frame:skip]):
        frame = universe.trajectory.frame
        fig = plt.figure(figsize = [10, 10])
        # obtain cluser sizes in dict
        cluster_sizes = {}
        for cluster in set(clusters[frame_idx]):
            cluster_sizes[cluster] = np.asarray(
                    np.nonzero(clusters[frame_idx] == cluster)).shape[1]
        # remove 0 cluster
        print(cluster_sizes)
        cluster_sizes.pop(0)
        cluster_counter = 0
        for cluster_key in cluster_sizes:
            if cluster_sizes[cluster_key] >= min_cluster_size:
                cluster_counter += 1
        fig.suptitle('Frame {} at {} ns, with {} clusters.'.format(
                frame, 
                frame*universe.trajectory.dt/1000, 
                cluster_counter))
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        ax.set_xlim3d(0, universe.dimensions[0])
        ax.set_ylim3d(0, universe.dimensions[1])
        ax.set_zlim3d(0, universe.dimensions[2])
        with open('cluster_sizes_over_time.dat', 'a') as f:
            f.write('{}\t{}\n'.format(frame, cluster_sizes))
        sorted_clusters = sorted(list(set(clusters[frame_idx])))
        for cluster in sorted_clusters:
            if cluster == 0:
                continue
            ax.scatter(universe.atoms.positions[clusters[frame_idx] == cluster][:,0][::reduce_points], universe.atoms.positions[clusters[frame_idx] == cluster][:,1][::reduce_points], universe.atoms.positions[clusters[frame_idx] == cluster][:,2][::reduce_points], alpha = 0.5, color = colors[(cluster-1)%len(colors)])
        #plt.show()
        fig.savefig('figs_sorted/leaflets_frame-{:09d}.png'.format(universe.trajectory.frame), dpi = 300)
        plt.close()
  

if __name__=='__main__':
    try:
        import clustering_input as inp
    except ModuleNotFoundError:
        print('There should be a file called clustering_input.py with needed settings. (An exmaple file should be made here)')
        sys.exit()
    start_frame = inp.start_frame
    stop_frame = inp.stop_frame
    skip = inp.skip

#    colors = get_spaced_colors(30)
#    colors = np.array(colors)
#    #TODO this has to be changed
#    colors = colors/255
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']


    print('Loading universe...')
    data = mda.Universe(inp.tpr, inp.xtc)
    print('Loading clusters...')
    clusters = np.load('clusters_ordered.npy')
    # this should be in atoms
    min_size = 150
    plot_clusters(data, clusters, skip, 1, colors, start_frame, stop_frame, min_cluster_size = min_size)
        
