# coding: utf-8
import numpy as np
import MDAnalysis as mda
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
import time
import itertools

def gen_explicit_matrix(atomgroup, resolution = 1, PBC = 'cubic'):
    """
    Takes an atomgroup and bins it as close to the resolution as possible.
    PBC is 'cubic' by default, but can be turned off. No other form of 
    PBC is currently supported.
    
    Returns
    (array) 3d boolean with True for occupied bins
    (dictionary) atom2voxel mapping
    """
    # scaling from ansgtrom to nm
    positions = atomgroup.positions/10
    # obtaining the matrix raw dimensions
    dimensions = atomgroup.dimensions[:3]/10
    
    # rounding the dimensions to the closest multiple of the resolution
    mod_dimensions = np.array(divmod(dimensions, resolution))
    mod_dimensions[1] = np.round(mod_dimensions[1]/resolution)
    mod_dimensions[0] = mod_dimensions[0]+mod_dimensions[1]
    
    # scaling the matrix to the binning
    scaling = mod_dimensions[0]/(dimensions/resolution)
    max_error = np.max(np.absolute(scaling-1))
    if max_error > 0.05:
        raise ValueError('A scaling artifact has occured of more than 5%, {}% deviation from the target resolution in frame {} was detected. You could consider increasing the resolution.'.format(np.abs(scaling-1)*100, atomgroup.universe.trajectory.frame))
    scaled_positions = ((positions * scaling) / resolution).astype(int)
  
    # fixing cubic PBC
    if PBC == 'cubic':
        scaled_positions = scaled_positions % mod_dimensions[0]
        
    # making an empty explicit matrix
    explicit_matrix = np.zeros((mod_dimensions[0].astype(int)),)
    
    # filling the explicit_matrix
    scaled_positions = scaled_positions.astype('int')
    explicit_matrix[scaled_positions[:,0], scaled_positions[:,1], scaled_positions[:,2]] = 1
    
    # generating the mapping dictionary
    voxel2atom = collections.defaultdict(list)
    for atom_index, scaled_position in enumerate(scaled_positions):
        x, y, z = scaled_position
        voxel2atom['x{}y{}z{}'.format(x, y, z)].append(atom_index)
    
    return explicit_matrix.astype(bool), voxel2atom

def blur_matrix(matrix, span = 1, PBC = 'cubic'):
    """
    Blurs a 3d boolean matrix by adding the values which are within span range.
    Default behaviour is using the inverse of the input matrix. This results in 
    an inner smearing useful for generating the inner contour. The opisite is 
    useful for generating the outer smearing for the outer contour. By default
    it assumes cubic periodic boundary conditions. PBC can currently not be 
    turned off. 
    """
    blurred_matrix = copy.copy(matrix)
    if PBC == 'cubic':
        for shift in range(span):
            for current_axis in range(3):
                blurred_matrix += np.roll(matrix, shift+1, current_axis)
                blurred_matrix += np.roll(matrix, -(shift+1), current_axis)
    return blurred_matrix.astype(bool)
   
def gen_contour(matrix, span = 1, inv = True):
    """
    Generatates the inner (inv = True), or outer (inv = False) contour of span 
    around the 3d boolean matrix.
    """
    if inv:
        blurred_matrix = blur_matrix(np.logical_not(matrix), span)
        return blurred_matrix ^ np.logical_not(matrix)
    if not inv:
        blurred_matrix = blur_matrix(matrix, span)
        return blurred_matrix ^ matrix

def find_neighbours(position, span):
    """
    Uses the position to generate a a box width size span around the position.
    """
    return itertools.product(
            range(position[0]-span, position[0]+span+1),
            range(position[1]-span, position[1]+span+1),
            range(position[2]-span, position[2]+span+1),
            )
    
# The 3d example of a very clear more set oriented neighbour clustering
def set_clustering(explicit_matrix, exclusion_mask = False, span = 1, verbose = False):
    """
    The 3d example of a very clear more set oriented neighbour voxel clustering.
    Generation of the set is relatively slow and could maybe be further optimized.
    
    Input should be a 3d boolean array.
    
    Returns a dictionary of clusters with a list of voxels per cluster.
    
    !!!TO BE DONE!!!
    Add an extra 3d bool input array for exclusion (useful for contaminated clusters such 
       as membrane leaflets with proteins).
    Add an extra input dictionary of clusters in the previous frame (useful for
       simple cluster identity over time). This can be used to draw from, and assign
       designated cluster based on the precious frame. A split will be indicated by an
       an overwrite in the dictionary for the cluster name if it has been split up in 
       the current frame.
       
       Example: (cluster, frame)
       # split (determenistic)
       (1,0)         --> [(11,1), (12,1)]
       # merge (determenistic, but not order independent)
       [(1,0),(2,0)] --> (1,1) or (2,1)
       # preserved identity
       (1,0)         --> (1,1)
       # lost to void
       (1,0)         --> cluster is no longer used from now on
       # new from void ($NEXT_AVAIL_CLUSTER is the highest ever used cluster plus one)
       (-,0)         --> ($NEXT_AVAIL_CLUSTER, 1)
    """
    # Obtaining a set for all occupied voxels
    # starting the timer for generating the set (verbose)
    start = time.time()
    # finding all hits (occupied voxels)
    if exclusion_mask is not False:
        positions = np.array(np.where(explicit_matrix ^ exclusion_mask == 1)).T
    else:
        positions = np.array(np.where(explicit_matrix == 1)).T
    if verbose:
        print('There are {} points to cluster.'.format(positions.shape[0]))
    # initiating the empty set for all hits
    positions_set = set()
    # adding each hit to the set as a tuple
    for position in positions:
        position = tuple(position)
        positions_set.add((position))
    # stop timer for making the hits set (verbose)
    stop = time.time()
    if verbose:
        print('It took {} to make the set.'.format(stop-start))

    # The clustering scales linear and millions of points can be achieved in the minute range.
    # beginning the timer for clustering (verbose)
    start = time.time()
    # the range of the neighbour search (1 is direct neighbour including the diagonal)
    span = 1
    # the starting cluster
    current_cluster = 1
    # output dictionary containing a list of hits per cluster
    clusters = {}
    # the to do queue
    queue = set()
    # getting the perdic dimensions
    dimensions = np.array(explicit_matrix.shape)
    # as long as there are hits
    while len(positions_set) > 0:
        # start the first point and remove from the hits
        current_position = positions_set.pop()
        # add self as first to current cluster
        clusters[current_cluster] = [current_position]
        # find all neighbours of self
        current_neighbours = list(itertools.product(range(int(current_position[0]-span), int(current_position[0]+span+1)),
                                                range(int(current_position[1]-span), int(current_position[1]+span+1)),
                                                range(int(current_position[2]-span), int(current_position[2]+span+1)),))
        # checking periodic boundary conditions
        for idx, neighbour in enumerate(current_neighbours):
            current_neighbours[idx] = (neighbour[0]%dimensions[0], neighbour[1]%dimensions[1], neighbour[2]%dimensions[2])
        # add all neighbours to the queue if they are in the hits
        queue =  positions_set.intersection(current_neighbours)
        while len(queue) > 0:
            # obtain current neighbour and remove from queue
            current_position = queue.pop()
            # also remove current neighbour from the hits
            positions_set.remove(current_position)
            # add self to current cluster
            clusters[current_cluster].append(current_position)
            # find all neighbours of self
            current_neighbours = list(itertools.product(range(int(current_position[0]-span), int(current_position[0]+span+1)),
                                                        range(int(current_position[1]-span), int(current_position[1]+span+1)),
                                                        range(int(current_position[2]-span), int(current_position[2]+span+1)),))
            # checking periodic boundary conditions
            for idx, neighbour in enumerate(current_neighbours):
                current_neighbours[idx] = (neighbour[0]%dimensions[0], neighbour[1]%dimensions[1], neighbour[2]%dimensions[2])
            # add all neighbours to the queue which are in the hits
            queue = queue.union(positions_set.intersection(current_neighbours))
        # move to next cluster
        current_cluster += 1
    # stopping the timer for clustering (verbose)
    stop = time.time()

    if verbose:
        print('It took {} to cluster {} clusters with a total of {} points'.format(stop-start, len(clusters), len(positions)))
    return clusters
    
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
    edge_color = (1,1,1,0.1)
    ax.voxels(array, edgecolor=edge_color, facecolor= color)
    plt.show()


if __name__=='__main__':
    data = mda.Universe('/home/bart/projects/clustering/test_files/4_adhesion/attached.gro')
    selection = data.select_atoms('resname DOPE DOTAP')
    resolution = 1

    start = time.time()
    explicit_matrix, voxel2atom = gen_explicit_matrix(selection, resolution = resolution)
    contour_matrix = gen_contour(explicit_matrix, 1, True)
    print('Making the contour took {}.\nGenerating output figures...'.format(time.time()-start))
    print('\nCLUSTERING 3, list based cubic boundary fix')
    clusters = set_clustering(contour_matrix, exclusion_mask = False, span = 1, verbose = True)
    plot_voxels(explicit_matrix)
    plot_voxels(contour_matrix)