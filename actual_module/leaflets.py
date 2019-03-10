# coding: utf-8
import MDAnalysis as mda
from . import clustering as clus
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import time
import sys
import pmda.custom

#multiprocess
#import pmda.custom

### Only for testing
def plot_voxels(array):
    """
    Plots a 3d boolean array.
    """
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
        atomgroup, exclusion_mask = False, resolution = 1, span = 1, inv = True
        ):
    """
    Clusters an mda.atomgroup based on its contour contacts in the active frame.
    
    Returns:
    --------
    A dictionary containing the voxels per cluster
    A dictionary containing the voxel2atoms conversion
    """
    # Generating the binary explicit matrix
    explicit_matrix, voxel2atoms = clus.gen_explicit_matrix(atomgroup, 
                                                            resolution)
    # calculating the contour mask
    contour_mask = clus.gen_contour(explicit_matrix, span, inv)
    # clustering the contours
    contour_clusters = clus.set_clustering(contour_mask, exclusion_mask)
    #plot_voxels(contour_mask)
    return contour_clusters, voxel2atoms, explicit_matrix, contour_mask

#@profile
def volume_clustering(
        atomgroup, exclusion_mask = False, resolution = 1
        ):
    """
    Clusters an mda.atomgroup based on its volume contacts in the active frame.
    
    Returns:
    --------
    A dictionary containing the voxels per cluster
    A dictionary containing the voxel2atoms conversion
    """
    explicit_matrix, voxel2atoms = clus.gen_explicit_matrix(atomgroup, 
                                                            resolution)
    volume_clusters = clus.set_clustering(explicit_matrix, exclusion_mask)
    return volume_clusters, voxel2atoms, explicit_matrix
    
def universe_clusters(clusters, mapping, selection):
    """
    Translates the voxels to atom indices in the selection.
    
    Currently clusters can steal from each other.
    """
    universe_masks = []
    for cluster in clusters:
        indexes = [
                mapping['x{}y{}z{}'.format(voxel[0], voxel[1], voxel[2])] 
                for voxel in clusters[cluster]
                ]
        indexes = np.concatenate(indexes).astype(int)
        assert np.unique(indexes).shape == indexes.shape, 'Indices should appear only once'
        universe_masks.append(selection[indexes])
    return universe_masks

def leaflet_clustering(
        universe, lipids_selection, tails_selection, 
        exclusions_selection = False, 
        resolution = 1
        ):
    """
    Clusters each lipid leaflet in the universe based on the 
    tails and full lipids of the given lipids in the universe. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position, this is usefull when the bilayer contains proteins. 
    All selections should be MDAnalysis atom_groups. The output is an 
    array with a cluster value for each atom. Clustering is atributed
    per residue.
    """   
    # creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection:
        # protein volume mask
        explicit_matrix_exclusions = clus.gen_explicit_matrix(
                exclusions_selection, resolution)[0]
        # protein contour (O) mask
        outward_contour_exclusions = clus.gen_contour(
                explicit_matrix_exclusions, inv = False
                )
        # protein volume+contour(O) mask
        exclusion_mask = explicit_matrix_exclusions+outward_contour_exclusions
        exclusion_mask[exclusion_mask > 1] = 1
    else:
        exclusion_mask = False
        
    # clustering the tail density for tail grouping
    current_selection = tails_selection
    current_clusters, current_mapping, volume_mask = volume_clustering(
            current_selection, exclusion_mask, resolution
            )
    #print('Current clusters: {}'.format(current_clusters.keys()))
    tails_universe_masks = universe_clusters(current_clusters, 
                                             current_mapping, 
                                             current_selection)
    #print('Tail universe masks: {}'.format(tails_universe_masks))
    tail_density_resid_groups = [
            tails_universe_mask.residues
            for tails_universe_mask in tails_universe_masks
            ]
    #print('Tail density resid groups: {}'.format(tail_density_resid_groups))
    # clustering the lipid contours per tail density group 
    list_lipid_contour_resid_groups = []
    for tail_density_resid_group in tail_density_resid_groups:
        current_selection = tail_density_resid_group.atoms.intersection(lipids_selection)
        current_clusters, current_mapping, cluster_volume, cluster_contour = contour_clustering(
                current_selection, exclusion_mask, resolution
                )
        # Expanding contour per contour cluster
        #  this is a still experimental part, the voxels are not only grown 
        #  inwards, thus the set check might be more expensive. We could 
        #  first make an explicit matrix for the whole lipds in the current
        #  lipid tail density group, but I am not sure if this is less overhead.
        #  A double hit voxel should not be silently overwritten, but not 
        #  assigned at all. This is also not happening for the first layer
        #  of the leaflets. This does have a performance impact in the range of
        #  10% tested for the 4_transfection lipoplex.
        for cluster in current_clusters:
            #print(current_clusters[cluster].shape)
            # we do not want to expand the void cluster
            if cluster == 0:
                continue
            contour_mask = np.zeros(cluster_contour.shape, dtype = bool)
            # creating contour cluster specific mask
            selection = np.array(current_clusters[cluster])
            contour_mask[selection[:,0], selection[:,1], selection[:,2]] = True
            # smearing the contour cluster outwards
            contour_mask = clus.gen_contour(contour_mask, inv = False)
            # adding the new oxels to the cluster
            extra_voxels = []
            extra_voxels = np.array(np.where(contour_mask == 1)).T
            current_clusters[cluster] = np.vstack((
                    current_clusters[cluster], extra_voxels)).astype(int)

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
    # generating the final array per frame indicating for each atom its cluster
    out_array = np.zeros((len(universe.atoms), ), dtype = 'uint32')
    for idx, leaflet_selection in enumerate(leaflet_selections):
        out_array[leaflet_selection.indices] = idx+1
    
    return out_array

def mf_leaflet_clustering(universe, lipids_selection, 
                          tails_selection, exclusions_selection = False, 
                          resolution = 1, skip = 1,  
                          return_selections = True, 
                          ):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails_selection and lipids_selection in the corresponding universe. 
    The output is an array containing the cluster per atom per frame. 
    Skip is used to skip frames for analysis.
    """
 
    start = time.time()
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
        
        # the actual clustering
        universe_mask = leaflet_clustering(universe, lipids_selection, 
                                            tails_selection, 
                                            exclusions_selection,
                                            resolution) 
        clusters.append(universe_mask)  
        
    print()
    clusters = np.array(clusters, dtype = 'uint32')
    #print('Clusters > 0: {}'.format(clusters[clusters > 0]))
    return clusters

def main():
    print('Reading trajectory...')
    # some local test data
    try:
        output_file = sys.argv[3] 
        plotting = sys.argv[4]
        skip = int(sys.argv[5])
        data = mda.Universe(sys.argv[1], sys.argv[2])
    except IndexError:
        print('Please specify a tpr, xtc, output, a 0/1 for plotting and a skip integer.')
        sys.exit()

    plotting = bool(int(plotting))
    # some basic lipid stuff, this will be moved to an input file
    # lipid selection tether
    #lipids = ['PIPI', 'POPA', 'PIPS', 'DUPE', 'CHOL', 'IPC', 'XNG1', 'DAPE', 'DBSM', 'PQPE', 'DXG3', 'OPC', 'XIPC', 'PQPS', 'PEPC', 'DPCE', 'PPC', 'DUPS', 'PAPS', 'PIPA', 'XAPS', 'PIDG', 'DPG3', 'POP1', 'UPC', 'DPG1', 'PAPC', 'BNSM', 'DXCE', 'POPE', 'PUPI', 'POP3', 'POPI', 'DAPC', 'DPSM', 'XNCE', 'PAPI', 'PAPA', 'APC', 'DOPE', 'DXG1', 'DAPS', 'PUPA', 'PODG', 'XOPE', 'XAPE', 'PNSM', 'PUPE', 'POSM', 'PNCE', 'XNSM', 'POP2', 'XPSM', 'PUPC', 'PUPS', 'DXSM', 'POPC', 'PGSM', 'XOPC', 'PADG', 'XHOL', 'PIPE', 'PUDG', 'POPS', 'DOPC', 'PAPE', 'XNG3', 'PIPC', 'PNG3', 'PNG1']
    # lipid selection lipoplex-membranes
    lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS', 'LYPC', 'LYPS']
    # the plasmamembrane lipids
    #lipids = ['PIPX', 'PEPC', 'PAPC', 'DAPC', 'POPE', 'PAPE', 'DAPE', 'PUPE', 'DPSM', 'DPG1', 'DXG1', 'PNG1',	'XNG1',	'DPG3',	'PNG3',	'XNG3',	'PIDG',	'CHOL',	'PIPC',	'DUPE',	'PGSM',	'PNSM',	'PIPS',	'PAPS',	'PIPI',	'POPA',	'POP1',	'PUDG',	'POPX',	'DOPE',	'PIPE',	'DBSM',	'DXG3',	'DPCE',	'DXCE',	'PPC',  'POPC',	'PQPE',	'POPS',	'PUPS',	'XNSM',	'PNCE',	'PADG',	'DOPC',	'PUPC',	'DXSM',	'POSM',	'BNSM',	'XNCE',	'APC',  'PODG',	'POPI',	'DUPS',	'IPC',  'UPC',  'PQPS',	'PAPI',	'PUPI',	'POP2',	'POP3',	'DAPS',	'PUPA',	'PIPA',	'PAPA',	'OPC']
    lipids_query = ' '.join(lipids)
    lipids_selection = data.select_atoms('resname {}'.format(lipids_query))
    #lipids_selection = lipids_selection.select_atoms('name C1A C1B C2A C2B C3A C3B C4A C4B C5A C5B D1A D1B D2A D2B D3A D3B D4A D4B D5A D5B GL1 GL2 or (resname CHOL and name C1 C2 R1 R2 R3 R4 R5 ROH)')
    #headgroups = ['CNO', 'NC3', 'NH3', 'PO4', 'TAP', 'GL1', 'GL2']   
    tails = ['C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A', 'D3B']
    tails_query = ' '.join(tails)
    tails_selection = data.select_atoms('name {} or (resname CHOL and name C1 C2)'.format(tails_query))
    #exclusions_selection = ['NA+']
    #exclusions_selection = ['BB1', 'BB2', 'SC1', 'SC2','SC3', 'SC4']
    #exclusions_query = ' '.join(exclusions_selection)
    #exclusions_selection = data.select_atoms('name {}'.format(exclusions_query))
    exclusions_selection = False
    resolution = 1
    print('Actual clustering...')
    #do_all = pmda.custom.analysis_class(leaflet_clustering)
    #clusters = do_all(data.trajectory, data, lipids_selection, tails_selection, exclusions_selection, resolution).run(step=skip, n_blocks=1)
    clusters = mf_leaflet_clustering(data, lipids_selection, tails_selection, 
                                     exclusions_selection, 
                                     resolution, skip)
    # writing the output at once
    np.save(output_file, clusters.astype('uint32'))                                 
    
    # some basic plotting
    if plotting:
        try:
            os.mkdir('figs')
        except FileExistsError:
            pass
        print('Some basic plotting...')
        universe = data
        reduce_points = 29
        for frame_idx, _ in enumerate(universe.trajectory[::skip]):
            frame = universe.trajectory.frame
            fig = plt.figure(figsize = [10, 10])
            fig.suptitle('Frame {} at {} ns, with {} clusters.'.format(
                    frame, 
                    frame*universe.trajectory.dt/1000, 
                    len(set(clusters[frame_idx]))-1)
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
            #turn on plotting hardcore
            #plt.show()
            plt.close()
    
    return clusters

if __name__ == '__main__':
    main()