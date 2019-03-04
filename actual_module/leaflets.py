# coding: utf-8
import MDAnalysis as mda
from . import clustering as clus
import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import pickle
import time
import sys

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
        universe, exclusion_mask = None, resolution = 1, density = 0.01, 
        inv_density = False
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
    contour_cluster_state_mask, contour_clusters = clus.clustering(
            contour_mask, exclusion_mask)
    #plot_voxels(contour_mask)
    return contour_clusters, voxel2atoms

#@profile
def volume_clustering(
        universe, exclusion_mask = None, resolution = 1, density = 0.01, 
        inv_density = False
        ):
    """
    Clusters a selection based on its volume.
    """
    explicit_matrix, voxel2atoms = clus.generate_explicit_matrix(
        universe, resolution = resolution, density = density, 
        inv_density = inv_density, verbose = False
        )
    #plot_voxels(explicit_matrix)
    volume_cluster_state_mask, volume_clusters = clus.clustering(
            explicit_matrix, exclusion_mask)
    return volume_clusters, voxel2atoms
    
def universe_clusters(clusters, mapping, selection):
    universe_masks = []
    for cluster in range(1, len(clusters)):
        indexes = [
                mapping['x{}y{}z{}'.format(voxel[0], voxel[1], voxel[2])] 
                for voxel in clusters[cluster]
                ]
        indexes = np.concatenate(indexes)
        assert np.unique(indexes).shape == indexes.shape, 'Indices should appear only once'
        universe_masks.append(selection[indexes])
    return universe_masks

# An attempt of allowing clustering of leaflets containing proteins.
def leaflet_clustering4(
        universe, lipid_resnames, tail_names, exclusion_names = None, 
        resolution = 1, density = 0.01, inv_density = False
        ):
    """
    Clusters each lipid leaflet in the universe based on the the 
    tails and full lipids of the given lipids in the universe. It treats 
    lipids flanking the proteins as a cluster breaker, preventing clustering
    past their position. Lipid resnames, tail names and protein names should 
    be a list of strings. The output is a list of of leaflet atomgroups.
    """   
    # handling selection input
    lipids_query = ' or '.join(
            ['resname {}'.format(lipid) for lipid in lipid_resnames]
            )
    tails_query = ' or '.join(
            ['name {}'.format(tail) for tail in tail_names]
            )
    lipids_selection = universe.select_atoms(lipids_query)
    tails_selection = lipids_selection.select_atoms(tails_query)
    
    # creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusion_names is not None:
        exclusions_query = ' or '.join(
                ['name {}'.format(name) for name in exclusion_names]
                )
        exclusions_selection = universe.select_atoms(exclusions_query)
        # protein volume mask
        explicit_matrix_exclusions = clus.generate_explicit_matrix(
                exclusions_selection, resolution = resolution, 
                density = density, inv_density = inv_density, verbose = False
                )[0]

        # protein contour (O) mask
        outward_contour_exclusions = clus.smear_3d_matrix(
                explicit_matrix_exclusions, inv=False
                )
        # protein volume+contour(O) mask
        exclusion_mask = explicit_matrix_exclusions+outward_contour_exclusions
        exclusion_mask[exclusion_mask > 1] = 1
    else:
        exclusion_mask = None
        
    # clustering the tail density for tail grouping
    current_selection = tails_selection
    current_clusters, current_mapping = volume_clustering(
            current_selection, exclusion_mask, resolution
            )
    tails_universe_masks = universe_clusters(current_clusters, 
                                             current_mapping, 
                                             current_selection)
    tail_density_resid_groups = [
            tails_universe_mask.residues
            for tails_universe_mask in tails_universe_masks
            ]
 
    # clustering the lipid contours per tail density group 
    list_lipid_contour_resid_groups = []
    for tail_density_resid_group in tail_density_resid_groups:
        current_selection = tail_density_resid_group.atoms
        current_clusters, current_mapping = contour_clustering(
                current_selection, exclusion_mask, resolution
                )
        lipids_universe_masks = universe_clusters(current_clusters, 
                                                  current_mapping, 
                                                  current_selection)
        lipid_contour_resid_groups = [
                lipids_universe_mask.residues
                for lipids_universe_mask in lipids_universe_masks
                ]
        list_lipid_contour_resid_groups += lipid_contour_resid_groups
    
    # combining the contour and the density for leaflet clustering
    leaflet_selections = []
    for lipid_contour_resid_group in list_lipid_contour_resid_groups:
        for tail_density_resid_group in tail_density_resid_groups:
            current_residues = lipid_contour_resid_group.intersection(
                tail_density_resid_group
            )
            if current_residues:
                leaflet_selections.append(current_residues.atoms)
    universe_masks = leaflet_selections
    
    return universe_masks

def leaflet_clustering5(
        lipids_selection, tails_selection, exclusions_selection = None, 
        resolution = 1, density = 0.01, inv_density = False
        ):
    """
    Clusters each lipid leaflet in the universe based on the the 
    tails and full lipids of the given lipids in the universe. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position. All selections should be MDAnalysis atom_groups. 
    The output is a list of of leaflet atomgroups.
    """   
    # creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection is not None:
        # protein volume mask
        explicit_matrix_exclusions = clus.generate_explicit_matrix(
                exclusions_selection, resolution = resolution, 
                density = density, inv_density = inv_density, verbose = False
                )[0]
        # protein contour (O) mask
        outward_contour_exclusions = clus.smear_3d_matrix(
                explicit_matrix_exclusions, inv=False
                )
        # protein volume+contour(O) mask
        exclusion_mask = explicit_matrix_exclusions+outward_contour_exclusions
        exclusion_mask[exclusion_mask > 1] = 1
    else:
        exclusion_mask = None
        
    # clustering the tail density for tail grouping
    current_selection = tails_selection
    current_clusters, current_mapping = volume_clustering(
            current_selection, exclusion_mask, resolution
            )
    tails_universe_masks = universe_clusters(current_clusters, 
                                             current_mapping, 
                                             current_selection)
    tail_density_resid_groups = [
            tails_universe_mask.residues
            for tails_universe_mask in tails_universe_masks
            ]
 
    # clustering the lipid contours per tail density group 
    list_lipid_contour_resid_groups = []
    for tail_density_resid_group in tail_density_resid_groups:
        current_selection = tail_density_resid_group.atoms.intersection(lipids_selection)
        current_clusters, current_mapping = contour_clustering(
                current_selection, exclusion_mask, resolution
                )
        lipids_universe_masks = universe_clusters(current_clusters, 
                                                  current_mapping, 
                                                  current_selection)
        lipid_contour_resid_groups = [
                lipids_universe_mask.residues
                for lipids_universe_mask in lipids_universe_masks
                ]
        list_lipid_contour_resid_groups += lipid_contour_resid_groups
    
    # combining the contour and the density for leaflet clustering
    leaflet_selections = []
    for lipid_contour_resid_group in list_lipid_contour_resid_groups:
        for tail_density_resid_group in tail_density_resid_groups:
            current_residues = lipid_contour_resid_group.intersection(
                tail_density_resid_group
            )
            if current_residues:
                leaflet_selections.append(current_residues.atoms)
    universe_masks = leaflet_selections
    
    return universe_masks

def mf_leaflet_clustering(universe, lipids_selection, 
                          tails_selection, exclusions_selection, resolution = 1, 
                          density = 0.01, inv_density = False, 
                          plotting = False,
                          skip = 1, reduce_points = 10, 
                          return_selections = True):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails_selection and lipids_selection in the corresponding universe. 
    The output is a super list where each element has a list of leaflet 
    atomgroups. Skip is used to skip frames for analysis and reduce_points 
    reduces the amount of data plotted in the output. Return_selections can be 
    turned off if one would only like to render the output images, this should 
    be combined with plotting, this prevents the memory from filling up with 
    the selections in the list for parsing large files. 
    (only for true XXXL data)
    """
    if plotting:
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
        message = 'Frame {}/{}. Leaflet clustering will take {:.01f} more minutes.\r'
        print(message.format(universe.trajectory.frame,
                             len(universe.trajectory)-1,
                             time_total-time_working,
                             ),
                             end = '')
        sys.stdout.flush()
        # clustering single frame leaflets
        # Clustrering 1.0
        #universe_masks = leaflet_clustering(universe, lipid_resnames, 
        #                                    tail_names, resolution = resolution)
        universe_masks = leaflet_clustering5(lipids_selection, 
                                             tails_selection, 
                                             exclusions_selection,
                                             resolution = resolution)
        
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
                ax.scatter(universe_mask.atoms.positions[:,0][::reduce_points], universe_mask.atoms.positions[:,1][::reduce_points], universe_mask.atoms.positions[:,2][::reduce_points], alpha = 0.05)
            fig.savefig('figs/leaflets_frame-{:09d}.png'.format(universe.trajectory.frame), dpi = 300)
            #turn on plotting hardcore
            #plt.show()
            plt.close()

    if plotting:
        print()
    if return_selections:
        return clusters
    else:
        return

# generating test data MDA
def main():
    # some local test data
    try:
        data = mda.Universe(sys.argv[1], sys.argv[2])
        #selection_path = sys.argv[3] 
        plotting = sys.argv[3]
        skip = int(sys.argv[4])
    except IndexError:
        print('Please specify a tpr, xtc, a 0/1 for plotting and a skip integer.')
        sys.exit()

    plotting = bool(int(plotting))
    # some basic lipid stuff, this will be moved to an input file
    # lipid selection tether
    #lipids = ['PIPI', 'POPA', 'PIPS', 'DUPE', 'CHOL', 'IPC', 'XNG1', 'DAPE', 'DBSM', 'PQPE', 'DXG3', 'OPC', 'XIPC', 'PQPS', 'PEPC', 'DPCE', 'PPC', 'DUPS', 'PAPS', 'PIPA', 'XAPS', 'PIDG', 'DPG3', 'POP1', 'UPC', 'DPG1', 'PAPC', 'BNSM', 'DXCE', 'POPE', 'PUPI', 'POP3', 'POPI', 'DAPC', 'DPSM', 'XNCE', 'PAPI', 'PAPA', 'APC', 'DOPE', 'DXG1', 'DAPS', 'PUPA', 'PODG', 'XOPE', 'XAPE', 'PNSM', 'PUPE', 'POSM', 'PNCE', 'XNSM', 'POP2', 'XPSM', 'PUPC', 'PUPS', 'DXSM', 'POPC', 'PGSM', 'XOPC', 'PADG', 'XHOL', 'PIPE', 'PUDG', 'POPS', 'DOPC', 'PAPE', 'XNG3', 'PIPC', 'PNG3', 'PNG1']
    # lipid selection lipoplex-membranes
    #lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS', 'LYPC', 'LYPS']
    # the plasmamembrane lipids
    lipids = ['PIPX', 'PEPC', 'PAPC', 'DAPC', 'POPE', 'PAPE', 'DAPE', 'PUPE', 'DPSM', 'DPG1', 'DXG1', 'PNG1',	'XNG1',	'DPG3',	'PNG3',	'XNG3',	'PIDG',	'CHOL',	'PIPC',	'DUPE',	'PGSM',	'PNSM',	'PIPS',	'PAPS',	'PIPI',	'POPA',	'POP1',	'PUDG',	'POPX',	'DOPE',	'PIPE',	'DBSM',	'DXG3',	'DPCE',	'DXCE',	'PPC',  'POPC',	'PQPE',	'POPS',	'PUPS',	'XNSM',	'PNCE',	'PADG',	'DOPC',	'PUPC',	'DXSM',	'POSM',	'BNSM',	'XNCE',	'APC',  'PODG',	'POPI',	'DUPS',	'IPC',  'UPC',  'PQPS',	'PAPI',	'PUPI',	'POP2',	'POP3',	'DAPS',	'PUPA',	'PIPA',	'PAPA',	'OPC']
    lipids_query = ' '.join(lipids)
    lipids_selection = data.select_atoms('resname {}'.format(lipids_query))
    #lipids_selection = lipids_selection.select_atoms('name C1A C1B C2A C2B C3A C3B C4A C4B C5A C5B D1A D1B D2A D2B D3A D3B D4A D4B D5A D5B GL1 GL2 or (resname CHOL and name C1 C2 R1 R2 R3 R4 R5 ROH)')
    #headgroups = ['CNO', 'NC3', 'NH3', 'PO4', 'TAP', 'GL1', 'GL2']   
    tails = ['C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A', 'D3B']
    tails_query = ' '.join(tails)
    tails_selection = data.select_atoms('name {} or (resname CHOL and name C1 C2)'.format(tails_query))
    #exclusions_selection = None
    exclusions = ['BB', 'SC1', 'SC2', 'SC3', 'SC4']
    exclusions_query = ' '.join(exclusions)
    exclusions_selection = data.select_atoms('name {}'.format(exclusions_query))
    resolution = 1
    clusters = mf_leaflet_clustering(data, lipids_selection, tails_selection, exclusions_selection, 
                                     plotting = plotting, resolution = resolution, 
                                     skip = skip, reduce_points = 1, 
                                     return_selections = True)
    return clusters

if __name__ == '__main__':
    main()