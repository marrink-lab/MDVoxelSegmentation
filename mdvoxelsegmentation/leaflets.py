# coding: utf-8
import os
import time
import sys
import numpy as np
import multiprocessing as mp
#from functools import partial
import MDAnalysis as mda
import matplotlib.pyplot as plt
from . import clustering as seg
from mpl_toolkits.mplot3d import Axes3D


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


# #@profile
# def contour_segmentation(
#         atomgroup, exclusions_mask=False, span=0, inv=True, args):
#     """
#     Segments an mda.atomgroup based on its contour contacts in the active
#     frame.
    
#     Returns:
#     --------
#     A dictionary containing the voxels per segment.
#     A dictionary containing the voxel2atoms conversion.
#     The explicit matrix.
#     The contour matrix.
#     """
#     # Generating the binary explicit matrix
#     explicit_matrix, voxel2atoms, nbox = seg.gen_explicit_matrix(
#         atomgroup, 
#         args.resolution,
#         args.hyper_resolution,
#         )
#     # calculating the contour mask
#     contour_mask = seg.contour(explicit_matrix, nbox, span, inv)
#     # Segmenting the contours
#     contour_segments = seg.set_clustering(contour_mask, nbox, exclusions_mask)
#     #plot_voxels(contour_mask)
#     return contour_segments, voxel2atoms, explicit_matrix, contour_mask


#@profile
def volume_segmentation(atomgroup, headgroups_mask, exclusions_mask, args):
    """
    Segments an mda.atomgroup based on its volume contacts in the active frame.
    The headgroups selection can be used to remove lipid tail densities wich 
    occupy the same voxel as headgroups (useful for bilayers in close 
    proximity). The exclusions_mask makes occupied voxels act as a stop for
    segmentation.
    
    Returns:
    --------
    A dictionary containing the voxels per segment
    A dictionary containing the voxel2atoms conversion
    The explicit matrix.
    """
    explicit_matrix, voxel2atoms, nbox = seg.gen_explicit_matrix(
        atomgroup, 
        args.resolution,
        args.hyper_resolution,
        )
    
    if headgroups_mask is not False:
        explicit_matrix[headgroups_mask] = False
        
    volume_segments = seg.set_clustering(explicit_matrix, nbox, exclusions_mask)
    return volume_segments, voxel2atoms, explicit_matrix

def connected_components_segmentation(selection_headgroups_atomgroup, 
                                      selection_linkers_atomgroup,
                                      exclusions, 
                                      args):
    """
    Returns a single frame segmentation array with basic connected components
    segmentaion.
    
    Segments all connected components. It treats exclusions as a local stop, 
    preventing segmentation past their position.
    
    Force segmentation can be used to segment all lipid residues which are not 
    yet in a segment by picking the segment the residue is surrounded
    by most within the specified range. Also takes a minimum segment
    size in particles. Segments smaller than the cutoff will not be returned.

    Parameters
    ----------
    selection_headgroups_atomgroup : MDAnalysis.atomgroup
        The atomgroup used for connected components analysis.
    exclusions : MDAnalysis.atomgroup
        The atomgroup used as barrier for the connected components analysis.
    args : namespace
        All the arguments from the input argparsing.

    Returns
    -------
    A single frame segmentation array.

    """    
    # Creating the exclusion mask for segmenting around the exclusions
    #  this will be use to set the exclusion (flanking) pixels to touched
    #  in the segmentation queue. Therefore they will act as a stop. But 
    #  can be segmented themselves.
    if selection_exclusions_atomgroup:
        # Exclusion volume mask.
        explicit_matrix_exclusions, _, nbox = seg.gen_explicit_matrix(
            selection_exclusions_atomgroup, 
            args.resolution, 
            args.hyper_resolution,
            )
        # Exclusion contour (O) mask.
        outward_contour_exclusions = seg.contour(
                explicit_matrix_exclusions, nbox, span=1, inv=False
                )
        # Exclusion volume+contour(O) mask.
        exclusions_mask = np.logical_or(explicit_matrix_exclusions,
                                        outward_contour_exclusions)
    else:
        exclusions_mask = False
        
    # Segmenting the tail density for tail grouping excluding the
    #  tail densities which are masked by headgroups.
    headgroups_segments, headgroups_mapping, headgroups_mask = volume_segmentation(
            selection_headgroups_atomgroup, 
            False,
            exclusions_mask, 
            args,
            )
    
    # Converting the voxel mask to selection atom indices.
    headgroups_atomgroups = seg.clusters2atomgroups(
            headgroups_segments,
            headgroups_mapping,
            selection_headgroups_atomgroup,
            )
    
    # Using the atom indices to obtain residues in selection.
    headgroups_residuegroups = [
            headgroups_atomgroup.residues
            for headgroups_atomgroup in headgroups_atomgroups
            ]
    
    # Writing the segments, skipping 0.
    out_array = np.zeros((len(selection_headgroups_atomgroup.universe.atoms)),
                         dtype = args.bit_size)
    for segment, segment_resid_group in enumerate(
            headgroups_residuegroups):
        # Satisfying minimum segment size.
        complete_selection = segment_resid_group.atoms
        # The minimum size is already satisfied and does not need to be
        #  checked again.
        if len(complete_selection) > args.minimum_size:
            # Using the atom indices to write the segment in the
            #   universe.atoms array.
            out_array[complete_selection.ix] = segment+1

    
    # Tries to segment non segmented lipids to the segments surrounding their
    #  headgroups. Segment 0 is excluded.    
    
    if args.force_segmentation:
        seg.iterative_force_clustering(
            selection_linkers_atomgroup, 
            int(args.force_segmentation*(2/3)), 
            out_array, 
            np.unique(out_array), 
            max_cutoff=args.force_segmentation, 
            max_stop=args.recursion_depth, 
            cutoff_increment=1, 
            verbose=False)
        
        if args.force_info:
            non_segmented_atoms = seg.non_clustered(
                selection_headgroups_atomgroup,
                out_array)
            print('{} Non segmented particles after forced '
                  'segmentation with a cutoff of {} Angstrom:'.format(
                      non_segmented_atoms, args.force_segmentation))

    return out_array
    

def leaflet_segmentation(
        selection_headgroups_atomgroup, 
        selection_linkers_atomgroup, 
        selection_tails_atomgroup, 
        selection_exclusions_atomgroup, 
        args
        ):
    """
    Returns a single frame segmentation array.
    
    Segments each lipid leaflet. It treats 
    exclusions as a local segment stop, preventing segmentation
    past their position, this is usefull when the bilayer contains proteins. 
    
    All selections should be an MDAnalysis atomgroups. The output is an 
    array with a segment value for each atom. Segmentation is atributed
    per residue. The matrix is stored as a 'uint32' array by default.
    
    Force segmentation can be used to segment all lipid residues which are not 
    yet in a segment by picking the segment the residue is surrounded
    by most within the specified range. Also takes a minimum segment
    size in particles. Segments smaller than the cutoff will not be returned.
    """
    # Generating the explicit matrix of all headgroups for masking the 
    #  lipid tail densities.
    headgroups_mask, headgroups_mapping, nbox = seg.gen_explicit_matrix(
        selection_headgroups_atomgroup, 
        args.resolution, 
        args.hyper_resolution,
        )
    
    ### Creating the exclusion mask for segmenting around the proteins
    #   this will be use to set the protein (flanking) pixels to touched
    #   in the segmentation queue. Therefore they will act as a stop. 
    if selection_exclusions_atomgroup:
        # Protein volume mask.
        explicit_matrix_exclusions, _, nbox = seg.gen_explicit_matrix(
            selection_exclusions_atomgroup, 
            args.resolution, 
            args.hyper_resolution,
            )
        # Protein contour (O) mask.
        outward_contour_exclusions = seg.contour(
                explicit_matrix_exclusions, nbox, span=1, inv=False
                )
        # Protein volume+contour(O) mask.
        exclusions_mask = np.logical_or(explicit_matrix_exclusions,
                                        outward_contour_exclusions)
    else:
        exclusions_mask = False

    # Segmenting the tail density for tail grouping excluding the
    #  tail densities which are masked by headgroups.
    tails_segments, tails_mapping, tails_mask = volume_segmentation(
            selection_tails_atomgroup, 
            headgroups_mask,
            exclusions_mask, 
            args,
            )
    if args.verbose:
        print('Plotting the headgroup masked lipid tail mask.')
        plot_voxels(tails_mask)
    
    # Making the tail explicit matrix without the masking of headgroups.
    all_tails_mask = seg.gen_explicit_matrix(
            selection_tails_atomgroup,
            args.resolution,
            args.hyper_resolution
            )[0]
    if args.verbose:
        print('Plotting the unmasked lipid tail mask.')
        plot_voxels(all_tails_mask)

    # Converting the voxel mask to selection atom indices.
    tails_atomgroups = seg.clusters2atomgroups(
            tails_segments,
            tails_mapping,
            selection_tails_atomgroup,
            )
    # Using the atom indices to obtain residues in selection.
    tails_residuegroups = [
            tails_atomgroup.residues
            for tails_atomgroup in tails_atomgroups
            ]
   
    # Segmenting the lipid contours per tail density group.
    list_leaflet_residgroups = []
    for tails_residuegroup in tails_residuegroups:
        # Generating the explicit matrix for the headgroups in current tails.
        local_headgroupsatomgroup = (tails_residuegroup.atoms &
                                      selection_headgroups_atomgroup)
        headgroups_mask, headgroups_mapping, nbox = seg.gen_explicit_matrix(
            local_headgroupsatomgroup, 
            args.resolution,
            args.hyper_resolution,
            )
        if args.verbose:
            print('Plotting the headgroups mask.')
            plot_voxels(headgroups_mask)
        
        # Using the unmasked tail densities to exclude headgroup locations.
        #  This will result in not segmenting all headgroups. Only covering 
        #  tails in the current selection are used as exclusion mask.
        headgroups_mask[all_tails_mask &
                         headgroups_mask] = False
        # Segmenting the masked headgroup densities.
        headgroups_segments = seg.set_clustering(
            headgroups_mask, 
            nbox, 
            exclusions_mask,
            )
        if args.verbose:
            print('\nCurrent leaflet segments: {}'.format(len(
                    headgroups_segments)))
            for segment in headgroups_segments:
                temp_mask = np.zeros(headgroups_mask.shape, dtype=bool)
                selection = np.array(headgroups_segments[segment])
                temp_mask[selection[:,0], selection[:,1],
                          selection[:,2]] = True
                plot_voxels(temp_mask)
        
        # Converting the voxel mask to selection atom indices with respect to 
        #  the universe.atoms.
        leaflets_atomgroups = seg.clusters2atomgroups(
                headgroups_segments,
                headgroups_mapping,
                local_headgroupsatomgroup,
                )
        # Converting the atom indices in selection to residues in selection.
        leaflets_residuegroups = [
                leaflets_atomgroup.residues
                for leaflets_atomgroup in leaflets_atomgroups
                ]
        # Adding the current resid groups to the list.
        list_leaflet_residgroups += list(leaflets_residuegroups)
    
    # Writing the lipid_contour_resid_groups per segment, skipping 0.
    out_array = np.zeros((len(selection_headgroups_atomgroup.universe.atoms)),
                         dtype = args.bit_size)
    for segment, leaflet_resid_group in enumerate(
            list_leaflet_residgroups):
        # Satisfying minimum segment size.
        leaflet_selection = leaflet_resid_group.atoms
        # The minimum size is already satisfied and does not need to be
        #  checked again.
        if len(leaflet_selection) > args.minimum_size:
            # Using the atom indices to write the segment in the
            #   universe.atoms array.
            out_array[leaflet_selection.ix] = segment+1

    
    # Tries to segment non segmented lipids to the segments surrounding their
    #  headgroups. Segment 0 is excluded.    
    
    if args.force_segmentation:
        seg.iterative_force_clustering(
            selection_linkers_atomgroup, 
            int(args.force_segmentation*(2/3)), 
            out_array, 
            np.unique(out_array), 
            max_cutoff=args.force_segmentation, 
            max_stop=args.recursion_depth, 
            cutoff_increment=1, 
            verbose=False)
        
        if args.force_info:
            non_segmented_atoms = seg.non_clustered(
                selection_headgroups_atomgroup,
                out_array)
            print('{} Non segmented particles after forced '
                  'segmentation with a cutoff of {} Angstrom:'.format(
                      non_segmented_atoms, args.force_segmentation))

    return out_array

    

def mf_leaflet_segmentation(universe,
                            headgroups_selection,
                            linkers_selection,
                            tails_selection, 
                            exclusions_selection,
                            start_frame, 
                            stop_frame,
                            args,
                            ):
    """
    MultiFrame Leaflet Segmentation
    
    Calls the leaflet_segmentation over a given range of frames. This is also
    where the loading bar is implemented. Returns the potential partial
    segments array (due to multithreading).
    
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
    segments = []
    # Iterating over each frame in trajectory using skip.
    for time_counter, _ in enumerate(universe.trajectory[
            start_frame:stop_frame:args.stride]):
        # LOADING BAR
        time_total = ((((time.time()-start)/(time_counter+1))*len(
                universe.trajectory[start_frame:stop_frame:args.stride])+1))/60
        time_working = (time.time()-start)/60
        message = 'Segmentation will take {:.01f} more minutes.\r'
        print(message.format(time_total-time_working), end = '')
        sys.stdout.flush()
        
        # The actual segmenting of individual frames.
        # Leaflet segmentation
        if headgroups_selection and tails_selection:
            segmentation = leaflet_segmentation(headgroups_selection,
                                                linkers_selection,
                                                tails_selection,
                                                exclusions_selection,
                                                args,
                                                )
        # Connected components segmentation.
        elif headgroups_selection and not tails_selection:
            segmentation = connected_components_segmentation(
                headgroups_selection, 
                linkers_selection,
                exclusions, 
                args):
        segments.append(segmentation)
        
    # This print is needed to get out of the same line as the loading bar of 
    #  the single frame leaflet segmentation.    
    #print() # Adds a newline to get out of the loading bar line.
    segments = np.array(segments, dtype = args.bit_size)
    return segments


def mf_leaflets_threaded(current_thread, args):
    """
    Returns the segmentation array. Takes arguments from the argparser and 
    spawns multiple threads for the leaflet segmentation.
    """
    # Generating the universe.
    #print('Reading trajectory...')
    universe = mda.Universe(args.reference, args.trajectory)

    # Importing selection queries from input file.
    # HEADGROUPS
    if args.headgroups_selection_query == 'False':
        headgroups_selection = False
    else:
        headgroups_selection = universe.select_atoms(
                args.headgroups_selection_query)
    
    # LINKERS
    if args.linkers_selection_query == 'False':
        linkers_selection = False
    else:
        linkers_selection = universe.select_atoms(args.linkers_selection_query)
    
    # TAILS
    if args.tails_selection_query == 'False':
        tails_selection = False
    else:
        tails_selection = universe.select_atoms(args.tails_selection_query)

    # EXCLUSIONS
    if args.exclusions_selection_query == 'False':
        exclusions_selection = False
    else:
        exclusions_selection = universe.select_atoms(
                args.exclusions_selection_query)
    
    # Checking which part of the trajectory should be treated.
    frames_per_process = int(args.frames // args.threads)
    start_frame = int(args.begin + (frames_per_process * current_thread))
    stop_frame = int(start_frame + frames_per_process)
    if current_thread == args.threads - 1:
        stop_frame = int(args.end)
        
    # Staring the segmentation.
    segments = mf_leaflet_segmentation(
        universe,
        headgroups_selection,
        linkers_selection,
        tails_selection,
        exclusions_selection, 
        start_frame, 
        stop_frame, 
        args)
    
    return segments
