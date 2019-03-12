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
    max_size = array.shape
    fig.suptitle('Plot dimensions: {}'.format(max_size))
    ax.set_xlim(0,max_size[0])
    ax.set_ylim(0,max_size[1])
    ax.set_zlim(0,max_size[2])
    color = (0.5,0.5,0.5,0.3)
    edge_color = (1,1,1,0.1)
    ax.voxels(array, edgecolor=edge_color, facecolor= color)
    plt.show()
    
#@profile
def contour_clustering(
        atomgroup, exclusion_mask = False, resolution = 1, span = 0, inv = True
        ):
    """
    Clusters an mda.atomgroup based on its contour contacts in the active frame.
    
    !!!! BUG I have a feeling that the inv False is not working as it should. (Bart)
    
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
        atomgroup, headgroups_selection = False,
        exclusion_mask = False, resolution = 1
        ):
    """
    Clusters an mda.atomgroup based on its volume contacts in the active frame.
    The headgroups selection can be used to remove lipid tail densities wich 
    occupy the same voxel as headgroups (useful for bilayers in close 
    proximity). The exclusion_mask makes occupied voxels act as a stop for
    clustering.
    
    Returns:
    --------
    A dictionary containing the voxels per cluster
    A dictionary containing the voxel2atoms conversion
    """
    explicit_matrix, voxel2atoms = clus.gen_explicit_matrix(atomgroup, 
                                                            resolution)
    if headgroups_selection is not False:
        explicit_matrix[headgroups_selection] = False
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

#@profile
def leaflet_clustering(
        universe, lipids_selection, tails_selection, headgroups_selection,
        exclusions_selection = False, 
        resolution = 1, expand = False
        ):
    """
    Clusters each lipid leaflet in the universe based on the 
    tails and full lipids of the given lipids in the universe. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position, this is usefull when the bilayer contains proteins. 
    All selections should be MDAnalysis atom_groups. The output is an 
    array with a cluster value for each atom. Clustering is atributed
    per residue. If expend is true, the outer two contours are used for
    the leaflet selection.
    """   
    # creating the headgroup exclusion mask for the tails
    if headgroups_selection:
        # protein volume mask
        headgroups_mask = clus.gen_explicit_matrix(
                headgroups_selection, resolution)[0]
    else:
        headgroups_mask = False
    
    # creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection:
        # protein volume mask
        explicit_matrix_exclusions = clus.gen_explicit_matrix(
                exclusions_selection, resolution)[0]
        # protein contour (O) mask
        outward_contour_exclusions = clus.gen_contour(
                explicit_matrix_exclusions, span = 1, inv = False
                )
        # protein volume+contour(O) mask
        exclusion_mask = np.logical_or(explicit_matrix_exclusions, 
                                       outward_contour_exclusions)
    else:
        exclusion_mask = False
        
    # clustering the tail density for tail grouping
    current_selection = tails_selection
    current_clusters, current_mapping, volume_mask = volume_clustering(
            current_selection, headgroups_mask, 
            exclusion_mask = False, resolution = resolution
            )
    # converting the voxel mask to selection atom indices
    tails_universe_masks = universe_clusters(current_clusters, 
                                             current_mapping, 
                                             current_selection)
    # converting the atom indices in selection to residues in selection
    tail_density_resid_groups = [
            tails_universe_mask.residues
            for tails_universe_mask in tails_universe_masks
            ]
    
    # clustering the lipid contours per tail density group 
    list_lipid_contour_resid_groups = []
    for tail_density_resid_group in tail_density_resid_groups:
        current_selection = tail_density_resid_group.atoms.intersection(lipids_selection)
        current_clusters, current_mapping, cluster_volume, cluster_contour = contour_clustering(
                current_selection, exclusion_mask, resolution
                )
        ## TEST PRINTS AND VOXEL PLOT
        #for cluster in current_clusters:
        #    temp_mask = np.zeros(volume_mask.shape, dtype=bool)
        #    selection = np.array(current_clusters[cluster])
        #    temp_mask[selection[:,0], selection[:,1], selection[:,2]] = True
        #    plot_voxels(temp_mask)
        #print('\nCurrent leaflet clusters: {}'.format(len(current_clusters)))
        #plot_voxels(cluster_contour)
        
        #! Expanding contour per contour cluster
        #    this is a still experimental part, the voxels are not only grown 
        #    inwards, thus the set check might be more expensive. We could 
        #    first make an explicit matrix for the whole lipds in the current
        #    lipid tail density group, but I am not sure if this is less overhead.
        #    A double hit voxel should not be silently overwritten, but not 
        #    assigned at all. This is also not happening for the first layer
        #    of the leaflets. This does have a performance impact in the range of
        #    10% tested for the 4_transfection lipoplex.
        # !!! BUG !!! NOT FUNCTIONAL OUTPUT IS NOT APPENDED !!!
        if expand:
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
        
        # converting the voxel mask to selection atom indices
        lipids_universe_masks = universe_clusters(current_clusters, 
                                                  current_mapping, 
                                                  current_selection)
        # converting the atom indices in selection to residues in selection
        lipid_contour_resid_groups = [
                lipids_universe_mask.residues
                for lipids_universe_mask in lipids_universe_masks
                ]
        # adding the current resid groups to the list
        #! does this also work when there is only one contour per tail density? 
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
    
    print(leaflet_selections)
    # generating the final array indicating for each atom its cluster
    out_array = np.zeros((len(universe.atoms), ))
    for idx, leaflet_selection in enumerate(leaflet_selections):
        out_array[leaflet_selection.indices] = idx+1
    
    return out_array

def leaflet_clustering(
        universe, lipids_selection, tails_selection, headgroups_selection,
        exclusions_selection = False, 
        resolution = 1, expand = False
        ):
    """
    Clusters each lipid leaflet in the universe based on the 
    tails and full lipids of the given lipids in the universe. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position, this is usefull when the bilayer contains proteins. 
    All selections should be MDAnalysis atom_groups. The output is an 
    array with a cluster value for each atom. Clustering is atributed
    per residue. If expend is true, the outer two contours are used for
    the leaflet selection.
    """   
    # creating the headgroup exclusion mask for the tails
    if headgroups_selection:
        # protein volume mask
        headgroups_mask, headgroups_mapping = clus.gen_explicit_matrix(
                headgroups_selection, resolution)
    else:
        headgroups_mask = False
    
    # creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection:
        # protein volume mask
        explicit_matrix_exclusions = clus.gen_explicit_matrix(
                exclusions_selection, resolution)[0]
        # protein contour (O) mask
        outward_contour_exclusions = clus.gen_contour(
                explicit_matrix_exclusions, span = 1, inv = False
                )
        # protein volume+contour(O) mask
        exclusion_mask = np.logical_or(explicit_matrix_exclusions, 
                                       outward_contour_exclusions)
    else:
        exclusion_mask = False
        
    # clustering the tail density for tail grouping
    tails_clusters, tails_mapping, tails_mask = volume_clustering(
            tails_selection, headgroups_mask, 
            exclusion_mask = False, resolution = resolution
            )
    # converting the voxel mask to selection atom indices
    tails_universe_masks = universe_clusters(tails_clusters, 
                                             tails_mapping, 
                                             tails_selection)
    # converting the atom indices in selection to residues in selection
    tail_density_resid_groups = [
            tails_universe_mask.residues
            for tails_universe_mask in tails_universe_masks
            ]
    
   
    # clustering the lipid contours per tail density group 
    list_lipid_contour_resid_groups = []
    for tail_density_resid_group in tail_density_resid_groups:
        #current_selection = tail_density_resid_group.atoms.intersection(headgroups_selection).select_atoms('name PO4 NC3 NH3 CNO GL1 GL2')
        current_selection = headgroups_selection.atoms.intersection(tail_density_resid_group.atoms)
        #current_selection = current_selection.select_atoms('name C4A C4B D4A D4B C3A C3B D3A D3B')
        #plot_voxels(headgroups_mask)
        #volume_mask = clus.gen_explicit_matrix(current_selection, resolution)[0]
        headgroups_mask, headgroups_mapping = clus.gen_explicit_matrix(current_selection, resolution)
        all_tails_mask = clus.gen_explicit_matrix(tails_selection, resolution)[0]
        headgroups_mask[all_tails_mask] = False
        #plot_voxels(headgroups_mask)
        #plot_voxels(volume_mask)
        #headgroup_mask[volume_mask] = False
        #plot_voxels(headgroup_mask)
        #plot_voxels(volume_mask)
        #plot_voxels(headgroups_mask)
        headgroups_clusters = clus.set_clustering(headgroups_mask)
        #print()
        #print('headgroup_clusters: {}'.format(len(headgroups_clusters)))
        #current_clusters, current_mapping, cluster_volume = volume_clustering(
        #        current_selection, False, exclusion_mask, resolution
        #        )
        ## TEST PRINTS AND VOXEL PLOT
        #for cluster in current_clusters:
        #    temp_mask = np.zeros(volume_mask.shape, dtype=bool)
        #    selection = np.array(current_clusters[cluster])
        #    temp_mask[selection[:,0], selection[:,1], selection[:,2]] = True
        #    plot_voxels(temp_mask)
        #print('\nCurrent leaflet clusters: {}'.format(len(current_clusters)))
        #plot_voxels(cluster_contour)
        
        #! Expanding contour per contour cluster
        #    this is a still experimental part, the voxels are not only grown 
        #    inwards, thus the set check might be more expensive. We could 
        #    first make an explicit matrix for the whole lipds in the current
        #    lipid tail density group, but I am not sure if this is less overhead.
        #    A double hit voxel should not be silently overwritten, but not 
        #    assigned at all. This is also not happening for the first layer
        #    of the leaflets. This does have a performance impact in the range of
        #    10% tested for the 4_transfection lipoplex.
             # !!! BUG !!! NOT FUNCTIONAL OUTPUT IS NOT APPENDED !!!
#        if expand:
#            for cluster in headgroups_clusters:
#                #print(current_clusters[cluster].shape)
#                # we do not want to expand the void cluster
#                if cluster == 0:
#                    continue
#                contour_mask = np.zeros(headgroups_mask.shape, dtype = bool)
#                # creating contour cluster specific mask
#                selection = np.array(headgroups_clusters[cluster])
#                contour_mask[selection[:,0], selection[:,1], selection[:,2]] = True
#                # smearing the contour cluster outwards
#                contour_mask = clus.gen_contour(contour_mask, span = 1, inv = False)
#                plot_voxels(contour_mask)
#                print()
#                print('Length cluster before appending: {}'.format(len(headgroups_clusters[cluster])))
#                headgroups_clusters[cluster] += np.array(np.where(contour_mask) == True).T
#                print('Length cluster after appending: {}'.format(len(headgroups_clusters[cluster])))
        
        # converting the voxel mask to selection atom indices
        lipids_universe_masks = universe_clusters(headgroups_clusters, 
                                                  headgroups_mapping, 
                                                  headgroups_selection)
        # converting the atom indices in selection to residues in selection
        lipid_contour_resid_groups = [
                lipids_universe_mask.residues
                for lipids_universe_mask in lipids_universe_masks
                ]
        # adding the current resid groups to the list
        #! does this also work when there is only one contour per tail density? 
        list_lipid_contour_resid_groups += lipid_contour_resid_groups
    
    print()
    print(len(list_lipid_contour_resid_groups))
    # combining the contour and the density for leaflet clustering
    leaflet_selections = []
    for lipid_contour_resid_group in list_lipid_contour_resid_groups:
        leaflet_selections.append(lipid_contour_resid_group.atoms)
    
    print(leaflet_selections)
    # generating the final array indicating for each atom its cluster
    out_array = np.zeros((len(universe.atoms), ))
    for idx, leaflet_selection in enumerate(leaflet_selections):
        print(len(leaflet_selection.intersection(leaflet_selections[0])))
        print(len(leaflet_selection.intersection(leaflet_selections[1])))
        print(len(leaflet_selection.intersection(leaflet_selections[2])))
        print(len(leaflet_selection.intersection(leaflet_selections[3])))
        out_array[leaflet_selection.indices] = idx+1
    
    return out_array


def mf_leaflet_clustering(universe, lipids_selection, 
                          tails_selection, headgroups_selection = False, 
                          exclusions_selection = False, 
                          resolution = 1, skip = 1, expand= False,   
                          return_selections = True, 
                          ):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails_selection and lipids_selection in the corresponding universe. 
    The output is an array containing the cluster per atom per frame. 
    Skip is used to skip frames for analysis. The cluster array is in a
    unsigned 32 bit format by default, allowing for 4294967296 unique clusters.
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
                                            tails_selection, headgroups_selection, 
                                            exclusions_selection,
                                            resolution, expand) 
        clusters.append(universe_mask)  
        
    print()
    clusters = np.array(clusters, dtype = 'uint32')
    #print('Clusters > 0: {}'.format(clusters[clusters > 0]))

    return clusters

def main():
    print('Reading trajectory...')
    try:
        import clustering_input as inp
    except IndexError:
        print('There should be a file called clustering_input with needed settings. (An exmaple file should be made here)')
        sys.exit()
    # generating the universe
    data = mda.Universe(inp.tpr, inp.xtc)

    # importing selection queries from input file
    lipids_selection = data.select_atoms(inp.lipids_selection_query)
    tails_selection = data.select_atoms(inp.tails_selection_query)
    if inp.headgroups_selection_query == False:
        headgroups_selection = False
    else:
        headgroups_selection = data.select_atoms(inp.headgroups_selection_query)
    if inp.exclusions_selection_query == False:
        exclusions_selection = False
    else:
        exclusions_selection = data.select_atoms(inp.exclusions_selection_query)
    
    # setting some other variables
    plotting = inp.plotting
    skip = inp.skip
    output_file = inp.output_file
    resolution = inp.resolution
    expand = inp.expand
    reduce_points = inp.reduce_points
    
    # staring the clustering
    print('Actual clustering...')
    start = time.time()
    clusters = mf_leaflet_clustering(data, lipids_selection, tails_selection, 
                                     headgroups_selection, exclusions_selection, 
                                     resolution, skip, expand)
    print('Clustering took: {}'.format(time.time()-start))
    #! writing the output at once this should become a per frame write and append!
    np.save(output_file, clusters.astype('uint32'))                                 

    # some basic plotting
    if plotting:
        try:
            os.mkdir('figs')
        except FileExistsError:
            pass
        print('Some basic plotting...')
        universe = data
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
            #turn on plotting hardcore
            #plt.show()
            plt.close()
    return clusters

if __name__ == '__main__':
    main()