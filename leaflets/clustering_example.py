#!python3
# coding: utf-8
import MDAnalysis as mda
import MDVoxelClustering.MDVoxelClustering as clus
import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import pickle
import time

### Only for testing
def plot_voxels(array):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca(projection='3d')
    max_size = np.array(array.shape).max()
    ax.set_xlim(0,max_size)
    ax.set_ylim(0,max_size)
    ax.set_zlim(0,max_size)    
    color = (0.5,0.5,0.5,0.3)
    edge_color = (1,1,1,0.3)
    ax.voxels(array, edgecolor=edge_color, facecolor= color)
    plt.show()

#@profile
def contour_clustering(
        universe, resolution = 1, density = 0.01, 
        inv_density = False, min_cluster_size = 5
        ):
    """
    Clusters a frame of an mda.Universe.selection object using their contour.
    
    Returns:
    --------
    A dictionary containing the voxels per cluster
    A dictionary containing the voxel2atoms conversion
    """
    # Generating the binary explicit matrix
    explicit_matrix, voxel2atoms = clus.generate_explicit_matrix(
        universe, resolution = resolution, density = density, 
        inv_density = inv_density, verbose = False
        )
    # calculating the contour mask
    contour_mask = clus.smear_3d_matrix(explicit_matrix)
    # clustering the contours
    contour_cluster_state_mask, contour_clusters = clus.clustering(contour_mask)
    return contour_clusters, voxel2atoms

def volume_clustering(
        universe, resolution = 1, density = 0.01, 
        inv_density = False, min_cluster_size = 5
        ):
    """
    Clusters a selection based on its volume.
    """
    explicit_matrix, voxel2atoms = clus.generate_explicit_matrix(
        universe, resolution = resolution, density = density, 
        inv_density = inv_density, verbose = False
        )
    volume_cluster_state_mask, volume_clusters = clus.clustering(explicit_matrix)
    return volume_clusters, voxel2atoms
    
def universe_clusters(clusters, mapping, selection):
    universe_masks = []
    for cluster in range(1, len(clusters)):
        indexes = [
                mapping['x{}y{}z{}'.format(voxel[0], voxel[1], voxel[2])] 
                for voxel in clusters[cluster]
                ]
        indexes = np.concatenate(indexes)
        universe_masks.append(selection[indexes])
    return universe_masks

def leaflet_clustering(
        universe, lipid_resnames, tail_names, 
        resolution = 1, density = 0.01, inv_density = False, 
        min_cluster_size = 5
        ):
    """
    Clusters each lipid leaflet in the universe based on the the 
    tails and full lipids of the given lipids in the universe. Names and 
    resname should be a list of strings. The output is a list of of leaflet
    universes.
    """   
    lipids_query = ' or '.join(
            ['resname {}'.format(lipid) for lipid in lipid_resnames]
            )
    # not being used atm
    #headgroups_query = ' or '.join(
    #        ['name {}'.format(headgroup) for headgroup in headgroup_names]
    #        )
    tails_query = ' or '.join(
            ['name {}'.format(tail) for tail in tail_names]
            )
    lipids_selection = universe.select_atoms(lipids_query)
    # not being used atm
    #headgroups_selection = lipids_selection.select_atoms(headgroups_query)
    
    tails_selection = lipids_selection.select_atoms(tails_query)

    # playing around with leaflet selector and it seems to work
    # clustering the lipid contours
    current_selection = lipids_selection
    current_clusters, current_mapping = contour_clustering(current_selection, resolution)
    lipids_universe_masks = universe_clusters(current_clusters, current_mapping, current_selection)
    lipid_contour_resid_groups = [
            set(lipids_universe_mask.resids)
            for lipids_universe_mask in lipids_universe_masks
            ]
    # clustering the tail density for tail grouping
    current_selection = tails_selection
    current_clusters, current_mapping = volume_clustering(current_selection, resolution)
    tails_universe_masks = universe_clusters(current_clusters, current_mapping, current_selection)
    tail_density_resid_groups = [
            set(tails_universe_mask.resids)
            for tails_universe_mask in tails_universe_masks
            ]
    
    # combining the contour and the density for leaflet clustering
    leaflet_selections = []
    for lipid_contour_resid_group in lipid_contour_resid_groups:
        for tail_density_resid_group in tail_density_resid_groups:
            current_resids = np.array(list(lipid_contour_resid_group.intersection(tail_density_resid_group)), dtype=int)
            if len(current_resids) >= 1:
                leaflet_universe = universe.atoms.residues[current_resids]        
                leaflet_selections.append(leaflet_universe)
    universe_masks = leaflet_selections
    
    return universe_masks


def mf_leaflet_clustering(universe, lipid_resnames, 
                               tail_names, resolution = 1, 
                               density = 0.01, inv_density = False, 
                               min_cluster_size = 5, plotting = False,
                               skip = 1, reduce_points = 10, 
                               return_selections = True):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails and complete lipids in the universe. The resnames and names for
    selection should be a list of strings. The output is a super list where 
    each element has a list of leaflet universes. Skip is used to skip frames
    for analysis and reduce_points reduces the amount of data plotted in the
    output. Return_selections can be turned off if one would only like to
    render the output images, this should be combined with plotting, this
    prevents the memory from filling up with the selections in the list for
    parsing large files.
    """
    try:
        os.mkdir('figs')
    except FileExistsError:
        pass
    
    start = time.time()
    if return_selections:
        clusters = []
    for _ in universe.trajectory[::skip]:
        # LOADING BAR
        
        time_total = ((((time.time()-start)/(universe.trajectory.frame+1))*len(universe.trajectory)+1))/60
        time_working = (time.time()-start)/60
        message = '\rCurrently working at frame {}/{}. The leaflet clustering will roughly take {:.02f} more minutes.'
        print(message.format(universe.trajectory.frame, 
                             len(universe.trajectory),
                             time_total-time_working,
                             ),
                             end = '')
        # clustering single frame leaflets
        universe_masks = leaflet_clustering(universe, lipid_resnames, 
                                            tail_names, resolution = resolution)
        if return_selections:
            clusters.append(universe_masks) 
        
        if plotting:
            fig = plt.figure(figsize = [10, 10])
            fig.suptitle('Frame {} at {} ns, with {} clusters.'.format(
                                universe.trajectory.frame, 
                                universe.trajectory.frame*universe.trajectory.dt/1000, 
                                len(universe_masks)) 
                        )
            ax = fig.add_subplot(111, projection='3d', aspect='equal')
            ax.set_xlim3d(0, universe.dimensions[0])
            ax.set_ylim3d(0, universe.dimensions[1])
            ax.set_zlim3d(0, universe.dimensions[2])
            for idx, universe_mask in enumerate(universe_masks):
                ax.scatter(universe_mask.atoms.positions[:,0][::reduce_points], universe_mask.atoms.positions[:,1][::reduce_points], universe_mask.atoms.positions[:,2][::reduce_points], alpha = 0.1)
            fig.savefig('figs/leaflets_frame-{:09d}.png'.format(universe.trajectory.frame), dpi = 300)
            plt.close()

    if plotting:
        print()
    if return_selections:
        return clusters
    else:
        return

# generating test data MDA
test = mda.Universe(
        '/home/bart/projects/clustering/MDVoxelClustering/test_files/4_adhesion/md.tpr', 
        '/home/bart/projects/clustering/MDVoxelClustering/test_files/4_adhesion/full_length.xtc',
        )

# some basic parsing
lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS', 'LYPC', 'LYPS']
headgroups = ['CNO', 'NC3', 'NH3', 'PO4', 'TAP', 'GL1', 'GL2']
tails = ['C3A', 'C3B', 'C4A', 'C4B']

clusters = mf_leaflet_clustering(test, lipids, tails, plotting = True, 
                                 skip = 5, reduce_points = 20, 
                                 return_selections = True)

#not working, but very useful to make it work
#with open('cluster_list.pkl', 'wb') as f:
#    pickle.dump(clusters, f)

## weird syntaxt, but flexible in a sense, it generates a rename selection
##  from a list.
#selection_resnames = ['DOPE', 'DOTAP']
##selection_headgroups = ['PO4', '']
##selection_resnames = ['DLPC', 'DLPS', 'LYPC', 'LYPS', 'DYPC', 'DYPS']
#selection_string = ' or '.join([' '.join(subquery) for subquery in list(
#        itertools.product(['resname'], selection_resnames))]
#        )
#lipoplex_lipids = u.select_atoms(selection_string)
## unit conversion from Angstrom to nm.
#test_data = lipoplex_lipids.positions/10
#print('{} particles to cluster.'.format(test_data.shape[0]))
#plotting = False
#
### set universe active frame
##u.trajectory[45]
#
#total_clusters = []
#for _ in u.trajectory[:]:
#    print('\rWorking at frame {}/{}'.format(u.trajectory.frame, 
#                                            len(u.trajectory)-1), end='')
#    contour_clusters, voxel2atoms = contour_clustering(test_data)
#    total_clusters.append(len(contour_clusters)-1)
#print()
#print(total_clusters)


#### transform the compressed matrix into an explicit matrix
#explicit_matrix, voxel2atoms = clus.generate_explicit_matrix(
#        test_data, resolution = 1, density = 0.01, inv_density = False
#        )
#if plotting:
#    print('Plotting the density mask:')
#    clus.plot_voxels(explicit_matrix)
#
#### clustering the densities
#
#density_cluster_state_mask, density_clusters = clus.clustering(explicit_matrix) 
## plotting the density clusters
#if plotting:
#    print('The density cluster(s):')
#    clus.plot_clusters(
#            density_cluster_state_mask, density_clusters, min_cluster_size = 5
#            )


#### calculating the contour mask
#contour_mask = clus.smear_3d_matrix(explicit_matrix)
#if plotting:
#    print('Plotting the contour mask:')
#    clus.plot_voxels(contour_mask)
#
#### clustering the contours
#contour_cluster_state_mask, contour_clusters = clus.clustering(contour_mask)
## plotting the contour clusters
#if plotting:
#    print('The contour cluster(s):')
#    clus.plot_clusters(
#            contour_cluster_state_mask, contour_clusters, 
#            min_cluster_size = min_cluster_size
#            )
#
#print('{} density cluster(s) and {} contour cluster(s) were detected.'.format(
#        len(density_clusters)-1, len(contour_clusters)-1)
#        )
