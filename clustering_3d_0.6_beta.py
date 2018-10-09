
# coding: utf-8

# In[101]:


##### An attempt of making data analysis easier and quicker for many MD people
### MAIN GOAL
#   To do analysis of MD trajectory on a more abstract level. Instead of having to talk about specific atoms
# and their corresponding position and residue members. We often would like to talk about properties on a
# much more abstract level, such as the diffusion of a certain aggregate of particles. Or to talk about the
# leaflets of lipid assemblies. Current methods often imply either geometrical arguments which can quickly
# become very expensive. Here we present a mixture of voxel and graph based selection methods to obtain 
# information about dynamic clusters/coarse particles at a hopefully realtime speads.

### CONTOURS
#   Contours are meant for analysis on the behaviour of clusters in a cheaper dimensionality. They should be
# enough to uniquely identify e.g. the number of particles present and the voxel space contour size (not 
# exactly the real size, but often close enough).

### VOLUMES
#   Volumes give you a more robust manner of selection which will for sure include all the particles you need
# for high resolution analysis. Combining layers of volumes and contours often can dramatically reduce the 
# degrees of freedom in your data set. Without loosing any particles on the way. Hopefully we will demonstrate 
# this for the automatic leaflet detection in a crowded membrane space (all forms of lipid cluster states, 
# such as: adhesed, semi-fused, fused and seperated).

#from pbcpy.base import pbcarray # needed for pbc clustering
import numpy as np
import copy
#import nglview as nv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import itertools
import collections
import datetime
import sys

### initial settings
# sets the resolution of the bins the amount in (nm/resolution)*3 is the minimum cluster distance at
#  which seperation is detected
#resolution = float(sys.argv[1])
resolution = 1
#the density is the ratio of the mean of the occupancy of occupied bins
#density = float(sys.argv[2])
# The density can be useful to make a distinction between edge paticles and bulk.  
#  The inverse tag allows you to select only the not dense regions and not 0.
density = 0.01
inv_density = False
#blur_density = 1 # the amount of neighbours is 27 at most
# warning setting the smearing from hard boundary to soft, results in good clustering,
#  but probably bad underlying particle restults, use with care!!! (we are going to use this for leaflet
#  detection)
#  show plots or not
plotting = True
min_cluster_size = 1

### Data, later this will be taken from MDAnalysis
#test_data = np.genfromtxt('attached_tails.dat')
#test_data = np.genfromtxt('ball.dat')
#test_data = np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0],[2,0,0],[2,1,0],
#                      [2,2,0],[0,2,0],[1,2,0],[0,0,1],[0,1,1],[1,1,1],[1,0,1],
#                      [2,0,1],[2,0,1],[2,1,1],[2,2,1],[0,0,2],[2,0,1],
#                      [0,1,2],[1,1,2],[1,0,2],[2,0,2],[2,0,2],[2,1,2],[2,2,2],
#                      [0,2,2],[1,2,2],[0,2,1],[1,2,1],
#                      [10,10,0], [10,11,0], [5,3,0],[3,5,0], [4,4,0],[4,5,0], [15,15,0]])
test_data = np.genfromtxt('D:\\Data\\eq_last_single_3d.dat')
print('{} particles to cluster.'.format(test_data.shape[0]))

#### This is where it all happens
def generate_explicit_matrix(array, resolution, density, 
                             specified_dim = False, inv_density = False, no_zero = True):
    """Takes a compressed 3d matrix and returns it as an explicit 3d matrix.
    The resolution is the relative bin size. A tuple of 3 can be used to speify the
    binning dimensions  in nm. It assumes your box has cubic PBC! Density can be used to 
    specify a minimum voxel density to be added to the output matrix. The inv_density can
    be set to true to specify a maximum density. The no_zero flag will even under the inv_density
    setting not return the elements containing 0 elements."""
    # we want to refer to the unmodified array for the dictionary entries later on
    array_original = copy.copy(array)
    # find the extremes to determine the final size of the explicit binned matrix
    x_max, y_max, z_max = np.max(array[:,0]), np.max(array[:,1]), np.max(array[:,2])
    if not specified_dim:
        limits = np.array([x_max, y_max, z_max])
    else:
        limits = np.array(specified_dim)
    # adapting the binning to the resolution
    limits = limits/resolution
    limits_ints = np.array(np.round(limits), dtype=int)
    # making the explicit matrix
    explicit_matrix = np.zeros(limits_ints)
    # converting the data points
    array_original = array
    array = copy.copy(array_original)
    array = array/resolution # convert data to bins
    # creating a dicitonary with the atoms inside and the xyz coodirdiantes as keys
    voxel2atoms = collections.defaultdict(list)
    # adding each poin to the explicit matrix also cubic PBC!!!
    for idx, point in enumerate(array):
        # warning this implements cubic PBC!!! A similar trick can be done for others
        x = round(point[0] % (limits_ints[0]-1)).astype(int)
        y = round(point[1] % (limits_ints[1]-1)).astype(int)
        z = round(point[2] % (limits_ints[2]-1)).astype(int)
        #print(x, type(x), y, type(y), z ,type(z))
        explicit_matrix[x, y, z] += 1
        # mapping atoms to voxels
        key = 'x{}y{}z{}'.format(x, y, z)
        voxel2atoms[key].append(array_original[idx])
    mean = explicit_matrix[explicit_matrix > 0].flatten().mean()
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
    occupancy_mask = np.array(np.logical_not(occupancy_mask),dtype=float)
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
    blurred_matrix[blurred_matrix < 0] = 0 # the amount of neighbours is 27 at most    
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

    
#def clustering_PBC_cubic(mask_cluster_state_mask, to_do_lis, pbc_voxels):
#    """Fixes cubic pbc by checking the boundary planes for particles and doing a very cheap neighbour search
#    at the other side of the PBC plane. Then it adds the found pbc partners to the queue.
#    This could probably be made more general to work for all PBC variants in 3d space."""
#    selected_plane = np.array(np.where(mask_cluster_state_mask[0] == 1)).T
#    for pbc_voxel in selected_plane:
#        to_do_list.append(border_voxel[0]-1 % , border_voxel[1], border_voxel[2])
#    return selected_plane
        

# getting the inner logic solid for clustering
def clustering(array, distance = 1, pbc = True):
    """Takes an index (xyz tuple) and uses the neighbour_mask (array) to search for neighbours in the
    edge_cluster_mask. It also sets the the touched flag in place in the edge_cluster_state_mask. For all found
    and processed neighbours
    
    ### The matrix need to have no edges in its outer boundary therefore we cheat and add on to all x y z
    ### we need to do this before
    """
    start = datetime.datetime.now()
    mask_cluster_state_mask, mask_indices = clustering_preparing(array)
    # the beginning of the outer cluster loop
    counter = 0
    #time = 0 
    max_time = len(mask_indices)
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
    clusters = list(set(mask_cluster_state_mask[1].flatten())) # a set is always returned low to high?
    finish = datetime.datetime.now()
    print('It took {} (days:hours:seconds:decimals) to do the clustering.'.format(finish-start))
    print('{} cluster(s) have been found:'.format(len(clusters)-1))
    return mask_cluster_state_mask, clusters

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
    
            
def plot_clusters(array, clusters, min_cluster_size = 5):
    """Creates a voxel plot for the clusters in the array never plots cluster 0 and only
    shows clusters of size equal or larger than the minimum."""
    edge_color = np.array((1,1,1,0.3), dtype=float)
    color = np.array((1,1,1,0.3), dtype=float)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca(projection='3d')
    max_size = np.array(array.shape).max()
    ax.set_xlim(0,max_size)
    ax.set_ylim(0,max_size)
    ax.set_zlim(0,max_size)
    counter = 0
    colors = []
    for cluster in clusters:
        # automagically skips cluster 0 :D
        plot_array = copy.copy(array[1])
        plot_array[plot_array != cluster] = 0
        plot_array[plot_array > 0] = cluster
        color[:3] = np.random.rand(3)
        if np.count_nonzero(plot_array.flatten()) >= min_cluster_size:
            colors.append(copy.copy(color))
            counter += 1
            ax.voxels(plot_array, edgecolor=edge_color, facecolors = color)
    print('{} cluster(s) have been found >= {} (min_cluster_size)'.format(counter, min_cluster_size))
    plt.show()
    counter = 0
    for cluster in clusters:
        plot_array = copy.copy(array[1])
        plot_array[plot_array != cluster] = 0
        plot_array[plot_array > 0] = cluster
        if np.count_nonzero(plot_array.flatten()) >= min_cluster_size:
            color = colors[counter]
            counter += 1
            print('Cluster {}'.format(counter))
            fig = plt.figure(figsize=(10, 10))
            ax = fig.gca(projection='3d')
            max_size = np.array(array.shape).max()
            ax.set_xlim(0,max_size)
            ax.set_ylim(0,max_size)
            ax.set_zlim(0,max_size)
            ax.voxels(plot_array, edgecolor=edge_color, facecolors = color)
            plt.show()
    
### transform the compressed matrix into an explicit matrix
explicit_matrix, voxel2atoms = generate_explicit_matrix(test_data, resolution, 
                                                        density, inv_density = inv_density)
if plotting:
    print('Plotting the density mask:')
    plot_voxels(explicit_matrix)
print('\nThough hard to see, the density of the lipoplex is connected to the bilayer through its pbc upper limit.\n')
    
### doing the clustering for the volumes
print('The density cluster(s):')
density_cluster_state_mask, density_clusters = clustering(explicit_matrix) 
# plotting the density clusters
plot_clusters(density_cluster_state_mask, density_clusters, min_cluster_size)
print('\nThe density of the lipoplex is connected to the bilayer through its pbc upper limit. '
      'This causes the volumes to be one cluster. If one would like to seperate such cases, selecting only '
      'the lipid tails will result in two seperate entities.\n')

### calculating the contour particles
contour_mask = smear_3d_matrix(explicit_matrix)
if plotting:
    print('Plotting the contour mask:')
    plot_voxels(contour_mask)
print('\nAs a result from the connected densities, the contour of the lipoplex is connected to the bilayer '
      'lower leafet through its pbc upper limit. '
      'This causes those contours to be one cluster. If one would like to seperate such cases, selecting only '
      'the lipid tails will result in two separate entities.\n')

### doing the clustering for the contours
print('The contour cluster(s):')
contour_cluster_state_mask, contour_clusters = clustering(contour_mask)
# plotting the contour clusters
plot_clusters(contour_cluster_state_mask, contour_clusters, min_cluster_size)
print('\nAs we can see, we can get a correct selecion for the inner channels of the lipopoplex. We are also '
      'capable to seperate the two leaflets. We also cluster over PBC, however, the phosphate densities '
      'of two non-fused approximate bilayers is not resolved correctly to allow for leaflet detection in '
      'such cases. We are currently working on a more advanced procedure to allow for such cases. '
      'Nevertheless, we believe that in many cases the above demonstrated procedure could help a lot! \n\n'
      'Basically, its an overhaul of MDAnalysis.leaflet. An implementation of the contour part on top '
      'of the MDA.leaflet algorythm could/should achieve the same result.')


# In[102]:


print(contour_clusters)


# In[103]:


test = np.array((np.where(contour_cluster_state_mask[1] == 2))).T


# In[104]:


selected_atoms = []
for i in test:
    for j in voxel2atoms['x{}y{}z{}'.format(i[0], i[1], i[2])]:
        selected_atoms.append(j)


# In[105]:


# the first 10 partciles in cluster 2
selected_atoms[:10]

