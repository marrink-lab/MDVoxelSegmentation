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
    
    
def plot_clusters(universe, clusters, skip = 100, reduce_points = 1, 
                  min_size = 0, start_frame = 0, stop_frame = None):
    """
    Makes a scatter plot for each cluster in each frame. The colour of the
    cluster is dependent on the integer of the cluster. Makes a 'figs' folder 
    in the active directory and writes succesively numbered *.png files. The
    recude points reduces the amount of points plotted this has a huge effect
    on the writing time needed. If a cluster is smaller than the min_size, it
    is not taken into account.
    """
    try:
        os.mkdir('figs')
    except FileExistsError:
        pass
    print('Some basic plotting...')
    if stop_frame == None:
        stop_frame = len(universe.trajectory)
    for frame_idx, _ in enumerate(
            universe.trajectory[start_frame:stop_frame:skip]):
        cluster_sizes = [(x,list(clusters[frame_idx]).count(x)) for x in 
                         set(clusters[frame_idx])]
        amount_clusters = len([ y for x, y in cluster_sizes 
                               if y >= min_size ])-1
        
        frame = universe.trajectory.frame
        fig = plt.figure(figsize = [10, 10])
        fig.suptitle('Frame {} at {} ns, with {} clusters.'.format(
                frame, 
                frame*universe.trajectory.dt/1000, 
                amount_clusters))
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        ax.set_xlim3d(0, universe.dimensions[0])
        ax.set_ylim3d(0, universe.dimensions[1])
        ax.set_zlim3d(0, universe.dimensions[2])
        #print('Current frame: {}'.format(frame_idx))
        for cluster, count in cluster_sizes:
            if cluster == 0:
                continue
            if count < min_size:
                continue
            #print(len(clusters[frame_idx][clusters[frame_idx] == cluster]))
            ax.scatter(universe.atoms.positions[
                    clusters[frame_idx] == cluster][:,0][::reduce_points], 
universe.atoms.positions[clusters[frame_idx] == cluster][:,1][::reduce_points], 
universe.atoms.positions[clusters[frame_idx] == cluster][:,2][::reduce_points], 
alpha = 0.5)
        fig.savefig('figs/leaflets_frame-{:09d}.png'.format(
                universe.trajectory.frame), dpi = 300)
        plt.close()
    return
  
          
#@profile
def contour_clustering(
        atomgroup, exclusion_mask = False, resolution = 1, span = 0, 
        inv = True):
    """
    Clusters an mda.atomgroup based on its contour contacts in the active
    frame.
    
    Returns:
    --------
    A dictionary containing the voxels per cluster.
    A dictionary containing the voxel2atoms conversion.
    The explicit matrix.
    The contour matrix.
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
    The explicit matrix.
    """
    explicit_matrix, voxel2atoms = clus.gen_explicit_matrix(atomgroup, 
                                                            resolution)
    if headgroups_selection is not False:
        explicit_matrix[headgroups_selection] = False
    volume_clusters = clus.set_clustering(explicit_matrix, exclusion_mask)
    return volume_clusters, voxel2atoms, explicit_matrix


def leaflet_clustering(
        tails_selection, headgroups_selection,
        exclusions_selection = False, 
        resolution = 1, bits = 'uint32', verbose = False,
        ):
    """
    Clusters each lipid leaflet in the the tails and headgroups. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position, this is usefull when the bilayer contains proteins. 
    All selections should be MDAnalysis atomgroups. The output is an 
    array with a cluster value for each atom. Clustering is atributed
    per residue. The matrix is stored as a 'uint32' array by default.
    """
    # Generating the explicit matix of all headgroups for masking the 
    #  lipid tail densities.
    headgroups_mask, headgroups_mapping = clus.gen_explicit_matrix(
            headgroups_selection, resolution)
    
    ### Creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection:
        # Protein volume mask.
        explicit_matrix_exclusions = clus.gen_explicit_matrix(
                exclusions_selection, resolution)[0]
        # Protein contour (O) mask.
        outward_contour_exclusions = clus.gen_contour(
                explicit_matrix_exclusions, span = 1, inv = False
                )
        # Protein volume+contour(O) mask.
        exclusion_mask = np.logical_or(explicit_matrix_exclusions, 
                                       outward_contour_exclusions)
    else:
        exclusion_mask = False
        
    # Clustering the tail density for tail grouping excluding the 
    #  tail densities which are masked by headgroups.
    tails_clusters, tails_mapping, tails_mask = volume_clustering(
            tails_selection, headgroups_mask, 
            exclusion_mask = False, resolution = resolution
            )
    if verbose:
        print('Plotting the headgroup masked lipid tail mask.')
        plot_voxels(tails_mask)
    # Making the tail explicit matrix without the masking of headgroups. 
    all_tails_mask = clus.gen_explicit_matrix(tails_selection, resolution)[0]
    if verbose:
        print('Currently plotting the unmasked lipid tail mask.')
        plot_voxels(all_tails_mask)
    
    # Converting the voxel mask to selection atom indices.
    tails_atomgroup_masks = clus.convert_clusters2atomgroups(tails_clusters, 
                                                             tails_mapping, 
                                                             tails_selection)
    # Using the atom indices to obtain residues in selection.
    tails_residuegroup_masks = [
            tails_atomgroup_mask.residues
            for tails_atomgroup_mask in tails_atomgroup_masks
            ]
   
    # Clustering the lipid contours per tail density group.
    list_leaflet_resid_groups = []
    for tails_residuegroup_mask in tails_residuegroup_masks:
        #TODO This is curerntly being tested to fix the indexing bug.
        current_selection = tails_residuegroup_mask.atoms
        # Generating the explicit matrix for the full lipids.
        headgroups_mask, headgroups_mapping = clus.gen_explicit_matrix(
                current_selection, resolution)
        if verbose:
            print('Currently plotting the headgroups mask.')
            plot_voxels(headgroups_mask)
        # Using the unmasked tail densities to exclude headgroup locations.
        #  This will result in not clustering all headgroups. 
        headgroups_mask[all_tails_mask] = False
        # Clustering the masked headgroup densities.
        headgroups_clusters = clus.set_clustering(headgroups_mask, 
                                                  exclusion_mask)

        ### TEST PRINTS AND VOXEL PLOT
        if verbose:
            print('\nCurrent leaflet clusters: {}'.format(len(
                    headgroups_clusters)))
            for cluster in headgroups_clusters:
                temp_mask = np.zeros(headgroups_mask.shape, dtype=bool)
                selection = np.array(headgroups_clusters[cluster])
                temp_mask[selection[:,0], selection[:,1], 
                          selection[:,2]] = True
                plot_voxels(temp_mask)
        
        # Converting the voxel mask to selection atom indices with respect to 
        #  the universe.atoms.
        leaflets_atomgroup_masks = clus.convert_clusters2atomgroups(
                headgroups_clusters, 
                headgroups_mapping, 
                current_selection
                )
        # Converting the atom indices in selection to residues in selection.
        leaflets_residuegroup_masks = [
                leaflets_atomgroup_mask.residues
                for leaflets_atomgroup_mask in leaflets_atomgroup_masks
                ]
        # Adding the current resid groups to the list.
        list_leaflet_resid_groups += list(leaflets_residuegroup_masks)
    print()
    
    # Writing the lipid_contour_resid_groups per cluster, skipping 0.
    out_array = np.zeros((len(headgroups_selection.universe.atoms)), 
                         dtype = bits)
    for cluster, lipid_contour_resid_group in enumerate(
            list_leaflet_resid_groups):
        leaflet_selection = lipid_contour_resid_group.atoms
        #TODO This is curerntly being tested to fix the indexing bug.
        # Using the atom indices to write the cluster in the universe.atoms 
        #   array.
        out_array[leaflet_selection.ix] = cluster+1
    return out_array
    

def mf_leaflet_clustering(universe, 
                          tails_selection, headgroups_selection = False, 
                          exclusions_selection = False, 
                          resolution = 1, skip = 1, bits = 'uint32',
                          verbose = False, start_frame = 0, stop_frame = None
                          ):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails_selection and headgroups_selection in the corresponding universe. 
    The output is an array containing the cluster per atom per frame. 
    Skip is used to skip frames for analysis. The cluster array is in a
    unsigned 32 bit format by default, allowing for 4294967296 unique clusters.
    """
    start = time.time()
    clusters = []
    if stop_frame == None:
        stop_frame = len(universe.trajectory)
    # Iterating over each frame in trajectory using skip.
    for time_counter, _ in enumerate(universe.trajectory[
            start_frame:stop_frame:skip]):
        # LOADING BAR
        time_total = ((((time.time()-start)/(time_counter+1))*len(
                universe.trajectory[start_frame:stop_frame:skip])+1))/60
        time_working = (time.time()-start)/60
        message = 'Frame {}/{}. Leaflet clustering will take {:.01f} more \
minutes.\r'
        print(message.format(universe.trajectory.frame,
                             stop_frame-1,
                             time_total-time_working,), 
                             end = '')
        sys.stdout.flush()
        
        # The actual clustering of individual frames.
        universe_mask = leaflet_clustering(tails_selection, 
                                           headgroups_selection, 
                                           exclusions_selection,
                                           resolution, bits, verbose) 
        clusters.append(universe_mask)  
    # This print is needed to get out of the same line as the loading bar of 
    #  the single frame leaflet clustering.    
    print() # Adds a newline to get out of the loading bar line.
    clusters = np.array(clusters, dtype = bits)
    return clusters


def main():
    # Parsing the input file.
    try:
        import clustering_input as inp
    except ModuleNotFoundError:
        print('There should be a file called clustering_input.py with needed\
settings. (An exmaple file should be made here)')
        sys.exit()
    # Generating the universe.
    print('Reading trajectory...')
    universe = mda.Universe(inp.tpr, inp.xtc)

    # Importing selection queries from input file.
    tails_selection = universe.select_atoms(inp.tails_selection_query)
    if inp.headgroups_selection_query == False:
        headgroups_selection = False
    else:
        headgroups_selection = universe.select_atoms(
                inp.headgroups_selection_query)
    if inp.exclusions_selection_query == False:
        exclusions_selection = False
    else:
        exclusions_selection = universe.select_atoms(
                inp.exclusions_selection_query)
    
    # Setting some other variables.
    plotting = inp.plotting
    skip = inp.skip
    output_file = inp.output_file
    resolution = inp.resolution
    reduce_points = inp.reduce_points
    verbose = inp.verbose
    start_frame = inp.start_frame
    stop_frame = inp.stop_frame
    #print('start {} stop {}'.format(start_frame, stop_frame))
    bits = 'uint32'
    
    # Staring the clustering.
    print('Actual clustering...')
    start = time.time()
    clusters = mf_leaflet_clustering(universe, tails_selection, 
                                     headgroups_selection, 
                                     exclusions_selection, 
                                     resolution, skip, bits, 
                                     verbose = verbose, 
                                     start_frame = start_frame, 
                                     stop_frame = stop_frame)
    print('Clustering took: {}'.format(time.time()-start))
    #TODO Writing the output at once this should become a per frame write/append!
    np.save(output_file, clusters.astype(bits))                                 

    # Some basic plotting.
    if plotting:
        plot_clusters(universe, clusters, skip, reduce_points, min_size = 150,
                      start_frame = start_frame, stop_frame = stop_frame)
    return clusters


if __name__ == '__main__':
    main()