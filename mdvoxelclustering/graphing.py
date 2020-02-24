#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import pickle


def find_total_amount_of_particles(frame, hidden_clusters=[0]):
    """
    Finds the total amount of non-zero particles.
    """
    # Finding the total amount of particles
    cluster_sizes_data[0]
    total_amount_of_particles = 0
    for cluster, count in cluster_sizes_data[0][frame]:
        if cluster in hidden_clusters:
            continue
        total_amount_of_particles += count
    return total_amount_of_particles


def main():
    # Reading in the clusters array
    print('Reading in cluster data.')
    clusters_array = np.load('clusters_ordered.npy')
    
    # Reading in clustering interactions pickle
    interactions_data = []
    with (open("cluster_mutations.pickle", "rb")) as openfile:
        while True:
            try:
                interactions_data.append(pickle.load(openfile))
            except EOFError:
                break
                
    # Reading in the cluster sizes pickle
    cluster_sizes_data = []
    with (open("visualization_data.pickle", "rb")) as openfile:
        while True:
            try:
                cluster_sizes_data.append(pickle.load(openfile))
            except EOFError:
                break
    
    
    print('Finding total amount of clusters.')
    # Finding highest cluster
    highest_cluster = -1
    for frame in clusters_array:
        unique_clusters = np.unique(frame)
        highest_local_cluster = max(unique_clusters)
        if highest_local_cluster > highest_cluster:
            highest_cluster = highest_local_cluster
    
    
    # Make the emtpy cluster counter dict
    cluster_size_dict = {}
    for cluster in range(highest_cluster+1):
        cluster_size_dict[cluster] = []
    
    
    number_of_frames = clusters_array.shape[0]
    # Fill up the dict
    for idx, frame in enumerate(clusters_array):
        print('\rProcessing frame {}/{} for instantiating the dictionary.'.format(idx, number_of_frames-1), end='')
        for cluster in cluster_size_dict:
            atom_count = len(np.where(frame == cluster)[0])
            cluster_size_dict[cluster].append(atom_count)
    print()
    
    
    # Smooth the cluster atom count trajectories
    N = 10
    for cluster in range(highest_cluster+1):
        print('\rProcessing cluster {}/{} for running average.'.format(cluster, highest_cluster), end='')
        cluster_size_dict[cluster] = np.array(cluster_size_dict[cluster])
        cluster_size_dict[cluster] = np.convolve(cluster_size_dict[cluster], np.ones((N,))/N, mode='same')
    print()
    
    
    # Obtain average cluster size and sort
    cluster_average_sizes = {}
    for cluster in cluster_size_dict:
        cluster_average_sizes[cluster] = np.average(cluster_size_dict[cluster])
        
    sorted_average_sizes_key_list = sorted(cluster_average_sizes, key=cluster_average_sizes.get, reverse=True)
    
    
    
    # Plot all clusters
    hidden_clusters = [0]
    max_clusters = 20
    counter = 0
    
    # Drawing the clusters
    fig, ax = plt.subplots()
    for idx, cluster in enumerate(sorted_average_sizes_key_list[:max_clusters]):
        if cluster in hidden_clusters:
            continue
        print('\rProcessing cluster {}/{} for plotting.'.format(idx+1, max_clusters), end='')
        ax.plot(cluster_size_dict[cluster], label='{}'.format(cluster))
    print()
    
    # Drawing the merge events between clusters weighted by size
    for frame in interactions_data[0]:
        for interaction in interactions_data[0][frame]:
            # We should be interested in this event based on its cluster id
            involved_cluster1 = interaction[0]
            if involved_cluster1 not in sorted_average_sizes_key_list[:max_clusters]:
                continue
            involved_cluster2 = interaction[1]
            if involved_cluster2 not in sorted_average_sizes_key_list[:max_clusters]:
                continue
            # We set up the line coordinates
            x1 = frame
            if involved_cluster1 in hidden_clusters:
                y1 = 0
            else:
                y1 = cluster_size_dict[involved_cluster1][frame]
            x2 = frame
            if involved_cluster2 in hidden_clusters:
                y2 = 0
            else:
                y2 = cluster_size_dict[involved_cluster2][frame]
            # We draw the connection line thickness related to the involved amount of particles
            total_amount_of_relevent_particles = find_total_amount_of_particles(frame, hidden_clusters)
            intensity = 0.9-0.9*(y1/total_amount_of_relevent_particles)
            ax.plot([x1, x2], [y1, y2], linestyle='-', linewidth=1, color=str(intensity))
    
    # Basic formatting and saving the graph
    fig.legend()
    ax.set_ylabel('atom count')
    ax.set_xlabel('frame number')
    plt.savefig('clusters_over_time.png', dpi=300)
    print('Plotting clusters complete.')

if __name__=='__main__':
    main()
