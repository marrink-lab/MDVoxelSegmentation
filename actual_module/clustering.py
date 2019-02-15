#!/usr/bin/env python
# coding: utf-8

#from pbcpy.base import pbcarray # needed for pbc clustering
#import MDAnalysis as mda
import numpy as np
import copy
#import nglview as nv
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
# allows for realtime plot rendering in notebook
#%matplotlib notebook
import itertools
import collections
import sys
import os

def read_config(config_file='MDclustering.inp', verbose = False):
    """
    A simple way to set default clustering settings.
    """
    file_is_present = os.path.isfile(config_file) 
    if file_is_present:
        exec(open(config_file).read())
        if verbose:
            with open(config_file) as f:
                document = f.readlines()
            print('The following settings will be used for the clustering.\n')
            for line in document:
                print(line.strip())
            print()
    else:
        print("Please create an input file (test_files/MDclustering.inp) with the following format:\n\nresolution  = 1\ndensity     = 0.01\ninv_density = False\nmin_cluster_size = 1\n")
        sys.exit()

#### This is where it all happens
#@profile
def generate_explicit_matrix(universe, resolution, density, specified_dim = False, 
                             inv_density = False, no_zero = True, verbose = False):
    """
    Takes a compressed 3d matrix and returns it as an explicit 3d matrix.
    The resolution is the relative bin size. A tuple of 3 can be used to 
    specify the binning dimensions  in nm. It assumes your box has cubic PBC! 
    Density can be used to specify a minimum voxel density to be added to the 
    output matrix. The inv_density can be set to true to specify a maximum 
    density. The no_zero flag will even under the inv_density setting not 
    return the elements containing 0 elements.
    """
    # protecting the original matrix
    # /10 for angstrom to nm conversion
    array = copy.copy(universe.positions/10)
    # find the extremes to determine the final size of the explicit binned matrix
    #x_max, y_max, z_max = np.max(array[:,0]), np.max(array[:,1]), np.max(array[:,2])
    if not specified_dim:
        limits = np.array([universe.dimensions[:3]/10])
    else:
        limits = np.array(specified_dim)
    # adapting the binning to the resolution
    limits = limits/resolution
    limits_ints = np.array(np.round(limits), dtype=int)
    limits_ints = limits_ints.flatten()
    # making the explicit matrix remove one for pbc
    explicit_matrix = np.zeros(limits_ints)
    # converting the data points
    array = array/resolution # convert data to bins
    # creating a dicitonary with the atoms inside and the xyz coodirdiantes as keys
    voxel2atoms = collections.defaultdict(list)
    # clipping the original matrix to the voxels
    # warning this implements cubic PBC!!! A similar trick can be done for others
    array = (array % (limits_ints.T)).astype(int)
    # adding each poin to the explicit matrix
    for idx, point in enumerate(array):
        x, y, z = point
        try:
            explicit_matrix[x, y, z] += 1
        except IndexError:
            print(limits_ints, point)
            return
        # mapping atoms to voxels
        key = 'x{}y{}z{}'.format(x, y, z)
        voxel2atoms[key].append(idx)
    mean = explicit_matrix[explicit_matrix > 0].flatten().mean()
    if verbose:
        print('The average bin density is {:.2f}'.format(mean))
    # clipping the matrix using the specified density or the inverse
    if inv_density:
        if no_zero:
            explicit_matrix[explicit_matrix == 0 ] = density*mean + 1
        explicit_matrix[explicit_matrix <= density*mean] = 1
        explicit_matrix[explicit_matrix > density*mean] = 0 
    else:
        explicit_matrix[explicit_matrix < density*mean] = 0 
        explicit_matrix[explicit_matrix >= density*mean] = 1
    return explicit_matrix, voxel2atoms

def smear_3d_matrix(array, pbc = True):
    """Takes an explicit array and smears it over the axes. This is like
    a running average in 3D. It returns the smeared array as a new array.
    By default the pbc is taken to be cubic, this can be turned off.
    
    !!! the smearing does not take pbc into account yet !!!"""
    shift = 1
    # making the matrix one bigger to prevent problems later on
    #dimensions = np.array(array.shape)+np.array((2, 2, 2))
    #array_empty = np.zeros(dimensions)
    #array_empty[1:-1,1:-1,1:-1] = array
    occupancy_mask = copy.copy(array)
    # inverting the matrix to get the inner boundaries
    occupancy_mask = np.array(np.logical_not(occupancy_mask),dtype=int)
    # smearing the matrix
    blurred_matrix = copy.copy(occupancy_mask)
    blurred_matrix[shift:] += occupancy_mask[:-shift]
    blurred_matrix[:-shift] += occupancy_mask[shift:]
    blurred_matrix[:,shift:] += occupancy_mask[:,:-shift]
    blurred_matrix[:,:-shift] += occupancy_mask[:,shift:]
    blurred_matrix[:,:,shift:] += occupancy_mask[:,:,:-shift]
    blurred_matrix[:,:,:-shift] += occupancy_mask[:,:,shift:]
    # smearing the pbc boundaries for cubic pbc
    if pbc:
        blurred_matrix[0] += occupancy_mask[-1]
        blurred_matrix[-1] += occupancy_mask[0]
        blurred_matrix[:,0] += occupancy_mask[:,-1]
        blurred_matrix[:,-1] += occupancy_mask[:,0]
        blurred_matrix[:,:,0] += occupancy_mask[:,:,-1]
        blurred_matrix[:,:,-1] += occupancy_mask[:,:,0]
    # clipping the matrix, tried density stuff here, but can't work.
    blurred_matrix[blurred_matrix >= 1] = 1   
    # obtaining the contours
    contour_mask = blurred_matrix - occupancy_mask
    return contour_mask

def clustering_preparing(array):
    mask_cluster_state_mask = np.array((array,
                                    np.zeros(array.shape),
                                    np.zeros(array.shape)))
    #mask_cluster_state_mask = pbcarray(mask_cluster_state_mask)
    temp_mask_indices = np.where(array == 1)
    mask_indices = np.array(list(zip(*temp_mask_indices)),dtype=int)
    
    return mask_cluster_state_mask, mask_indices

def cubic_selection(dimensions, position, distance = 1):
    """Creates a selection mask which can handle cubic pbc."""
    position = np.array(position)
    selection = np.array((position-distance, position+distance+1)).T
    #print('unmodified selection', selection)
    selection_list = []
    # create the non continuous axes for the 3 cases (too low,high, or normal)
    for axes, single_selection in enumerate(selection):
        #print('axes {} selection unmodified'.format(axes), single_selection)
        if single_selection[0] < 0:
            single_selection = (single_selection[1], dimensions[axes]+single_selection[0])
            temp_ones = np.ones(dimensions[axes])
            temp_ones[single_selection[0]:single_selection[1]] = 0
            #print(temp_ones)
            axes_selection = np.array(np.where(temp_ones)).flatten()
            #print('axes {} selection modified 1'.format(axes), single_selection)
            #print(axes_selection)
        elif single_selection[1] >= dimensions[axes]:
            single_selection = (single_selection[1] % (dimensions[axes]), single_selection[0])
            temp_ones = np.ones(dimensions[axes])
            temp_ones[single_selection[0]:single_selection[1]] = 0
            axes_selection = np.array(np.where(temp_ones)).flatten()
            #print('axes {} selection modified 2'.format(axes), single_selection)
            #print(axes_selection)
        else:
            axes_selection = np.arange(single_selection[0], single_selection[1]).flatten()
            #print('axes {} selection modified 3'.format(axes), single_selection)
            #print(axes_selection)
        selection_list.append(axes_selection)
    # create the coordinates for the neighbours in the non continous selection
    selection_mask = np.zeros(dimensions)
    for element in list(itertools.product(selection_list[0], selection_list[1], selection_list[2])):
        selection_mask[element] = 1
    return np.array(np.where(selection_mask == 1)).T

def clustering_inner_loop(idx, mask_cluster_state_mask, counter, to_do_list, distance, pbc = True):
    """A pretty smart clustering procedure."""
    # this automagically never flies out of bounds :D
    # this is also wher I have to implement PBC lets start with cbic
    #  the goal is to make the mask select the other side edges for negative numbers
    min_distance = idx - distance
    max_distance = idx + distance + 1
    simple_indices = np.logical_and(np.all(min_distance >= 0),
                            np.all(max_distance <= mask_cluster_state_mask.shape[1:]))
    # either no pbc, or the indexes are not an issue without pbc treatment
    if not pbc or simple_indices:
        neighbour_hits = mask_cluster_state_mask[:3, min_distance[0]:max_distance[0], 
                                                     min_distance[1]:max_distance[1], 
                                                     min_distance[2]:max_distance[2]]
        # This is where some cool masks are made for detection of neighbours
        neighbour_coordinates = neighbour_hits[0].astype(bool) # where?
        neighbour_clusters = neighbour_hits[1] # ref to cluster positions
        neighbour_clusters[neighbour_coordinates] = counter # set clusters for found neighbours
        neighbour_states = neighbour_hits[2].astype(bool) # check state of neighbour
            
        # if you are a neighbour and you have not been placed in the queue
        new_coordinates = np.where(neighbour_coordinates & ~neighbour_states)
        new_coordinates = np.array(new_coordinates).T
        try:
            # transform local to world coordinates
            new_coordinates = new_coordinates + idx - (1,1,1)
        except ValueError:
            return
        # add the global coordinates of the next iterations
        to_do_list.extend(new_coordinates)
        # setting all the neighbours to hit so they will never 
        #  be done again
        neighbour_hits[2, :] = 1
    elif pbc:
        neighbour_coordinates = cubic_selection(mask_cluster_state_mask.shape[1:], idx, distance)
        for neighbour_coordinate in neighbour_coordinates:
            # find occupied voxels this is either 0 or 1
            neighbour_hit = mask_cluster_state_mask[0,
                                                     neighbour_coordinate[0],
                                                     neighbour_coordinate[1],
                                                     neighbour_coordinate[2]]
            if neighbour_hit == 1:
                # set cluster 
                mask_cluster_state_mask[1,
                                         neighbour_coordinate[0],
                                         neighbour_coordinate[1],
                                         neighbour_coordinate[2]] = counter
                # append to queue if not touched before
                if mask_cluster_state_mask[2,
                                         neighbour_coordinate[0],
                                         neighbour_coordinate[1],
                                         neighbour_coordinate[2]] == 0:
                    to_do_list.append(neighbour_coordinate)
                # set touched state
                mask_cluster_state_mask[2,
                                         neighbour_coordinate[0],
                                         neighbour_coordinate[1],
                                         neighbour_coordinate[2]] = 1
    # this is probably not what the user had in mind
    else:
        print(idx, ' could not be processed with the current pbc settings.')

# getting the inner logic solid for clustering
def clustering(array, distance = 1, pbc = True):
    """Takes an index (xyz tuple) and uses the neighbour_mask (array) to search for neighbours in the
    edge_cluster_mask. It also sets the the touched flag in place in the edge_cluster_state_mask. For all found
    and processed neighbours
    
    ### The matrix need to have no edges in its outer boundary therefore we cheat and add on to all x y z
    ### we need to do this before
    """
    mask_cluster_state_mask, mask_indices = clustering_preparing(array)
    # the beginning of the outer cluster loop
    counter = 0
    #time = 0 
    for idx in mask_indices:
        # genreates a dequeue, which is like a list but cheaper tot pop at front 
        to_do_list = collections.deque()
        self_state = mask_cluster_state_mask[2, idx[0], idx[1], idx[2]]
        if self_state == 0:
            counter += 1
            # execute cluster function and start queue for cluster
            clustering_inner_loop(idx, mask_cluster_state_mask, counter, to_do_list, distance, pbc)
        while len(to_do_list) > 0: # exhaust all queue members for cluster
            idx = np.array(to_do_list.popleft())
            clustering_inner_loop(idx, mask_cluster_state_mask, counter, to_do_list, distance, pbc)
        
    clusters = set(mask_cluster_state_mask[1].flatten()) # a set is always returned low to high?
    cluster_dict ={}
    for x in range(len(clusters)):
        cluster_dict[x] = np.array(np.where(mask_cluster_state_mask[1,:,:,:] ==x)).T #Create a Dictionary to store voxels by cluster
    return mask_cluster_state_mask, cluster_dict


#  =============================================================================
#  def plot_voxels(array):
#      fig = plt.figure(figsize=(10, 10))
#      ax = fig.gca(projection='3d')
#      max_size = np.array(array.shape).max()
#      ax.set_xlim(0,max_size)
#      ax.set_ylim(0,max_size)
#      ax.set_zlim(0,max_size)    
#      color = (0.5,0.5,0.5,0.3)
#      edge_color = (1,1,1,0.3)
#      ax.voxels(array, edgecolor=edge_color, facecolor= color)
#      plt.show()
#  
#  def plot_clusters(array, clusters, min_cluster_size = 5):
#      """Creates a voxel plot for the clusters in the array never plots cluster 0 and only
#      shows clusters of size equal or larger than the minimum."""
#      edge_color = np.array((1,1,1,0.3), dtype=float)
#      color = np.array((1,1,1,0.3), dtype=float)
#  
#      fig = plt.figure(figsize=(10, 10))
#      ax = fig.gca(projection='3d')
#      max_size = np.array(array.shape).max()
#      ax.set_xlim(0,max_size)
#      ax.set_ylim(0,max_size)
#      ax.set_zlim(0,max_size)
#      counter = 0
#      colors = []
#      for cluster in clusters:
#          # automagically skips cluster 0 :D
#          plot_array = copy.copy(array[1])
#          plot_array[plot_array != cluster] = 0
#          plot_array[plot_array > 0] = cluster
#          color[:3] = np.random.rand(3)
#          if np.count_nonzero(plot_array.flatten()) >= min_cluster_size:
#              colors.append(copy.copy(color))
#              counter += 1
#              ax.voxels(plot_array, edgecolor=edge_color, facecolors = color)
#      print('{} cluster(s) have been found >= {} (min_cluster_size)'.format(counter, min_cluster_size))
#      plt.show()
#      counter = 0
#      for cluster in clusters:
#          plot_array = copy.copy(array[1])
#          plot_array[plot_array != cluster] = 0
#          plot_array[plot_array > 0] = cluster
#          if np.count_nonzero(plot_array.flatten()) >= min_cluster_size:
#              color = colors[counter]
#              counter += 1
#              print('Cluster {}'.format(counter))
#              fig = plt.figure(figsize=(10, 10))
#              ax = fig.gca(projection='3d')
#              max_size = np.array(array.shape).max()
#              ax.set_xlim(0,max_size)
#              ax.set_ylim(0,max_size)
#              ax.set_zlim(0,max_size)
#              ax.voxels(plot_array, edgecolor=edge_color, facecolors = color)
#              plt.show()
#  =============================================================================
