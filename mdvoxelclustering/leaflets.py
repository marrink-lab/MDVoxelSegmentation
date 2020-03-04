# coding: utf-8
import multiprocessing as mp
from functools import partial
import MDAnalysis as mda
from . import clustering as clus
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import time
import sys

mda.core.periodic = True

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
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim3d(0, universe.dimensions[0])
        ax.set_ylim3d(0, universe.dimensions[1])
        ax.set_zlim3d(0, universe.dimensions[2])
        print('Current frame: {}\r'.format(frame_idx), end='')
        for cluster, count in cluster_sizes:
            if cluster == 0:
                continue
            if count < min_size:
                continue
            #print(len(clusters[frame_idx][clusters[frame_idx] == cluster]))
            ax.scatter(
universe.atoms.positions[clusters[frame_idx] == cluster][:,0][::reduce_points], 
universe.atoms.positions[clusters[frame_idx] == cluster][:,1][::reduce_points], 
universe.atoms.positions[clusters[frame_idx] == cluster][:,2][::reduce_points], 
alpha = 0.5,
)
        fig.savefig('figs/leaflets_frame-{:09d}.png'.format(
                universe.trajectory.frame), dpi = 300)
        plt.close()
    print()
    return
  
          
#@profile
def contour_clustering(
        atomgroup, exclusion_mask=False, resolution=1, span=0, inv=True, 
        frames=0, hyper_res = False,
        ):
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
    explicit_matrix, voxel2atoms, nbox = clus.gen_explicit_matrix_multiframe(
        atomgroup, 
        resolution,
        frames=frames,
        hyper_res = hyper_res
    )
    # calculating the contour mask
    contour_mask = clus.contour(explicit_matrix, nbox, span, inv)
    # clustering the contours
    contour_clusters = clus.set_clustering(contour_mask, nbox, exclusion_mask)
    #plot_voxels(contour_mask)
    return contour_clusters, voxel2atoms, explicit_matrix, contour_mask


#@profile
def volume_clustering(
        atomgroup, headgroups_selection=False,
        exclusion_mask=False, resolution=1, frames=0, hyper_res = False,
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
    explicit_matrix, voxel2atoms, nbox = clus.gen_explicit_matrix_multiframe(
        atomgroup, 
        resolution,
        frames = frames,
        hyper_res = hyper_res,
    )
    if headgroups_selection is not False:
        explicit_matrix[headgroups_selection] = False
    volume_clusters = clus.set_clustering(explicit_matrix, nbox, exclusion_mask)
    return volume_clusters, voxel2atoms, explicit_matrix


def leaflet_clustering(
        selection_tails_atomgroup, selection_headgroups_atomgroup,
        exclusions_selection = False,
        resolution=1, bits='uint32', verbose=False, force=False, 
        force_cutoff=20, frames=1, force_info=True, hyper_res = False,
        min_cluster_size = 0
        ):
    """
    Clusters each lipid leaflet in the the tails and headgroups. It treats 
    exclusions as a local cluster stop, preventing clustering
    past their position, this is usefull when the bilayer contains proteins. 
    All selections should be MDAnalysis atomgroups. The output is an 
    array with a cluster value for each atom. Clustering is atributed
    per residue. The matrix is stored as a 'uint32' array by default.
    Force can be used to cluster all lipid residues which are not yet in a 
    cluster by picking the cluster the 0 bead from the residue is surrounded
    by most within the force_cutoff (Angstromg). Also takes a minimum cluster
    size in atoms. Clusters smaller than the cutoff will not be returned.
    """
    test = 1
    #TODO REMOVE PRINT
    if test:
        print('\n{}'.format(test))
        test += 1
    # Generating the explicit matix of all headgroups for masking the 
    #  lipid tail densities.
    headgroups_mask, headgroups_mapping, nbox = clus.gen_explicit_matrix_multiframe(
        selection_headgroups_atomgroup, resolution, frames=frames,
        hyper_res = hyper_res
    )
    
    #TODO REMOVE PRINT
    if test:
        print('{}'.format(test))
        test += 1
    ### Creating the exclusion mask for clustering around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the clustering queue. Therefore they will act as a stop. 
    if exclusions_selection:
        # Protein volume mask.
        explicit_matrix_exclusions, _, nbox = clus.gen_explicit_matrix_multiframe(
            exclusions_selection, resolution, frames=frames, 
            hyper_res = hyper_res
        )
        # Protein contour (O) mask.
        outward_contour_exclusions = clus.contour(
                explicit_matrix_exclusions, nbox, span=1, inv=False
                )
        # Protein volume+contour(O) mask.
        exclusion_mask = np.logical_or(explicit_matrix_exclusions,
                                       outward_contour_exclusions)
    else:
        exclusion_mask = False

    #TODO REMOVE PRINT
    if test:
        print('{}'.format(test))
        test += 1
    # Clustering the tail density for tail grouping excluding the
    #  tail densities which are masked by headgroups.
    tails_clusters, tails_mapping, tails_mask = volume_clustering(
            selection_tails_atomgroup, headgroups_mask,
            exclusion_mask=False, resolution=resolution, frames=frames
            )
    if verbose:
        print('Plotting the headgroup masked lipid tail mask.')
        plot_voxels(tails_mask)
    # Making the tail explicit matrix without the masking of headgroups.
    all_tails_mask = clus.gen_explicit_matrix_multiframe(
            selection_tails_atomgroup,
            resolution, frames=frames,
            hyper_res = hyper_res
            )[0]
    if verbose:
        print('Currently plotting the unmasked lipid tail mask.')
        plot_voxels(all_tails_mask)

    # Converting the voxel mask to selection atom indices.
    tails_atomgroups = clus.convert_clusters2atomgroups(
            tails_clusters,
            tails_mapping,
            selection_tails_atomgroup,
            frames,
            hyper_res,
            )
    # Using the atom indices to obtain residues in selection.
    tails_residuegroups = [
            tails_atomgroup.residues
            for tails_atomgroup in tails_atomgroups
            ]
   
    # Clustering the lipid contours per tail density group.
    list_leaflet_residgroups = []
    
    #TODO REMOVE PRINT
    if test:
        print('{}'.format(test))
        test += 1
        test_loop = 1
    
    for tails_residuegroup in tails_residuegroups:
        
        #TODO REMOVE PRINT
        if test:
            print('{}:{}'.format(test, test_loop))
            test_loop += 1
            
        # Generating the explicit matrix for the headgroups in current tails.
        local_headgroupsatomgroup = (tails_residuegroup.atoms &
                                      selection_headgroups_atomgroup)
        headgroups_mask, headgroups_mapping, nbox = clus.gen_explicit_matrix_multiframe(
            local_headgroupsatomgroup, resolution, frames=frames,
            hyper_res = hyper_res,
        )
        
        if verbose:
            print('Currently plotting the headgroups mask.')
            plot_voxels(headgroups_mask)
        
        # Using the unmasked tail densities to exclude headgroup locations.
        #  This will result in not clustering all headgroups. Only covering 
        #  tails in the current selection are used as exclusion mask.
        headgroups_mask[all_tails_mask &
                         headgroups_mask] = False
        # Clustering the masked headgroup densities.
        headgroups_clusters = clus.set_clustering(
            headgroups_mask, nbox, exclusion_mask
        )

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
        leaflets_atomgroups = clus.convert_clusters2atomgroups(
                headgroups_clusters,
                headgroups_mapping,
                local_headgroupsatomgroup,
                frames,
                hyper_res,
                )
        # Converting the atom indices in selection to residues in selection.
        leaflets_residuegroups = [
                leaflets_atomgroup.residues
                for leaflets_atomgroup in leaflets_atomgroups
                ]
        # Adding the current resid groups to the list.
        list_leaflet_residgroups += list(leaflets_residuegroups)
    
    # Writing the lipid_contour_resid_groups per cluster, skipping 0.
    out_array = np.zeros((len(selection_headgroups_atomgroup.universe.atoms)),
                         dtype = bits)
    for cluster, leaflet_resid_group in enumerate(
            list_leaflet_residgroups):
        # Satisfying minimum clustersize.
        leaflet_selection = leaflet_resid_group.atoms
        # The minimum size is already satisfied and does not need to be
        #  checked again.
        if len(leaflet_selection) > min_cluster_size:
            # Using the atom indices to write the cluster in the
            #   universe.atoms array.
            out_array[leaflet_selection.ix] = cluster+1

    
    # Tries to clusters non clustered lipids to the clusters surrounding their
    #  headgroups. Cluster 0 is excluded.    
    
    #TODO REMOVE THESE TESTS!
    if test:
        print('{}, now its recurion time'.format(test))
        test += 1
    
    #TODO Remove the hard codedness here.
    if force:
        clus.iterative_force_clustering(
            selection_headgroups_atomgroup.universe.select_atoms(
                '(name PO4 GL1 GL2 '
                'C1A C1B AM1 AM2 GM1 GM2 COO COOH) or '
                '(resname CHOL and name ROH) or '
                '(resname PAPI PIPI POP1 POP2 POP3 POPI PUPI and '
                'name C1 C2 C3 P1 P2 P3) or (name BB)'), 
            int(force_cutoff*(2/3)), out_array, 
            np.unique(out_array), max_cutoff=force_cutoff, max_stop=20, 
            cutoff_increment=1, verbose=False)
        
        if force_info:
            non_clustered_atoms = clus.non_clustered(
                selection_headgroups_atomgroup,
                out_array)
            print('{} Non clustered to be clustered atoms after forced '
                  'clustering with a cutoff of {} Angstrom:'.format(
                      non_clustered_atoms, force_cutoff))
    
    if test:
        print('{}, recursion done'.format(test))
        test += 1

    return out_array

    

def mf_leaflet_clustering(universe,
                          tails_selection, headgroups_selection=False,
                          exclusions_selection=False,
                          resolution=1, skip=1, bits='uint32',
                          verbose=False, start_frame=0, stop_frame=None,
                          force=True, force_cutoff=20, frames=1, 
                          force_info=True, hyper_res=False, 
                          min_cluster_size = 0):
    """
    MultiFrame Leaflet Clustering
    
    Clusters each lipid leaflet in the universe.trajectory based on the the 
    tails_selection and headgroups_selection in the corresponding universe. 
    The output is an array containing the cluster per atom per frame. 
    Skip is used to skip frames for analysis. The cluster array is in a
    unsigned 32 bit format by default, allowing for 4294967296 unique clusters.
    
    TODO I should not chunk the data in consecutive  chunks, its much better
    to perform the multiprocess at start+thread::num_threads. This allows
    for a much better performance in the parrelelization. This is due to the 
    fact there depending on the configuration in the box, the algorythm might
    have a harder time figuring it out. This in combination that there is a 
    time correlation in the data, this would mean that there are certain 
    consequtive regions which might be harder then others. Therefore it is
    wise to split the data in a skip manner, to divide correlated frames to
    different threads. There a much better total load balance is achieved.
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
                                           resolution, bits, verbose, force, 
                                           force_cutoff, frames=frames,
                                           force_info=force_info,
                                           hyper_res = hyper_res, 
                                           min_cluster_size = min_cluster_size 
                                           )
        clusters.append(universe_mask)
    # This print is needed to get out of the same line as the loading bar of 
    #  the single frame leaflet clustering.    
    print() # Adds a newline to get out of the loading bar line.
    clusters = np.array(clusters, dtype = bits)
    return clusters


def leaflet_clustering_threaded(
        current_frame,
        selection_tails_atomgroup, 
        selection_headgroups_atomgroup,
        exclusions_selection = False,
        resolution=1, bits='uint32', verbose=False, force=False, 
        force_cutoff=20, frames=1, force_info=True, hyper_res = False,
        min_cluster_size = 0
        ):
    """
    Threaded leaflet clutsering.
    """
    selection_tails_atomgroup.universe.trajectory[current_frame]
    clusters = leaflet_clustering(
            selection_tails_atomgroup,
            selection_headgroups_atomgroup,
            exclusions_selection,
            resolution, bits, verbose, force, 
            force_cutoff, frames,
            force_info, hyper_res = hyper_res,
            min_cluster_size = min_cluster_size
            )
    # The actual clustering of individual frames.
    clusters = np.asarray(clusters, dtype = bits)

    return clusters


def mf_leaflets_threaded(current_thread):
    # Parsing the input file.
    try:
        import clustering_input as inp
    except ModuleNotFoundError:
        print('There should be a file called clustering_input.py with needed\
settings. (An exmaple file should be made here)')
        sys.exit()
    # Generating the universe.
    #print('Reading trajectory...')
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
    resolution = inp.resolution
    hyper_res = inp.hyper_res
    min_cluster_size = inp.min_cluster_size
    force = inp.force
    force_cutoff = inp.force_cutoff
    verbose = inp.verbose
    frames = inp.frames
    force_info = inp.force_info
    start_frame = inp.start_frame
    stop_frame = inp.stop_frame
    skip = inp.skip
    threads = inp.threads
    bits = 'uint32'
    
    # multi processing bit
    if stop_frame == None:
        stop_frame = len(universe.trajectory)
    total_frames = stop_frame - start_frame
    frames_per_process = total_frames//threads
    start_frame = start_frame + (frames_per_process * current_thread)
    stop_frame = start_frame + frames_per_process
    if current_thread == threads - 1:
        stop_frame = inp.stop_frame
    # Staring the clustering.
    #print('Actual clustering...')
    clusters = mf_leaflet_clustering(
        universe,
        tails_selection, headgroups_selection,
        exclusions_selection,
        resolution, skip, bits,
        verbose, start_frame, stop_frame,
        force, force_cutoff, frames, 
        force_info, hyper_res = hyper_res,
        min_cluster_size = min_cluster_size)
    
    return clusters


def main_threaded():
    try:
        import clustering_input as inp
    except ModuleNotFoundError:
        print('There should be a file called clustering_input.py with needed\
settings. (An exmaple file should be made here)')
        sys.exit()
    # Generating the universe.
    print('Reading trajectory...')
    universe = mda.Universe(inp.tpr, inp.xtc)

    # Setting some other variables.
    plotting = inp.plotting
    skip = inp.skip
    reduce_points = inp.reduce_points
    start_frame = inp.start_frame
    stop_frame = inp.stop_frame
    bits = 'uint32'
    output_file = inp.output_file
    threads = inp.threads

    # Staring the clustering.
    print('Actual clustering...')
    start = time.time()
    
    pool = mp.Pool(threads)
    clusters = pool.map(mf_leaflets_threaded, range(threads))
    clusters = np.asarray(clusters)
    clusters = np.concatenate(clusters, axis=0)
    pool.close()
    pool.join()
    
    print(clusters, len(clusters))
    
    print('Clustering took: {}'.format(time.time()-start))
    #TODO Writing the output at once this should become a per frame 
    #  write/append!
    np.save(output_file, clusters.astype(bits))

    # Some basic plotting.
    if plotting:
        plot_clusters(universe, clusters, skip, reduce_points, min_size = 150,
                      start_frame = start_frame, stop_frame = stop_frame)

def main():
    # Parsing the input file.
    try:
        import clustering_input as inp
    except ModuleNotFoundError:
        print('There should be a file called clustering_input.py with needed '
              'settings. (An exmaple file should be made here)')
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
    hyper_res = inp.hyper_res
    min_cluster_size = inp.min_cluster_size
    force = inp.force
    force_cutoff = inp.force_cutoff
    reduce_points = inp.reduce_points
    verbose = inp.verbose
    start_frame = inp.start_frame
    stop_frame = inp.stop_frame
    frames = inp.frames
    force_info = inp.force_info
    bits = 'uint32'

    
    # Staring the clustering.
    print('Actual clustering...')
    start = time.time()
    clusters = mf_leaflet_clustering(
            universe, tails_selection,
            headgroups_selection,
            exclusions_selection,
            resolution, skip, bits,
            verbose=verbose,
            start_frame=start_frame,
            stop_frame=stop_frame, force=force, 
            force_cutoff=force_cutoff, frames=frames,
            force_info=force_info, hyper_res = hyper_res, 
            min_cluster_size = min_cluster_size,
            )
    print('Clustering took: {}'.format(time.time()-start))
    #TODO Writing the output at once this should become a per frame write/append!
    np.save(output_file, clusters.astype(bits))

    # Some basic plotting.
    if plotting:
        plot_clusters(universe, clusters, skip, reduce_points, min_size = 150,
                      start_frame = start_frame, stop_frame = stop_frame)


if __name__ == '__main__':
    main()
