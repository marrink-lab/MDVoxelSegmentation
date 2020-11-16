# coding: utf-8
import copy
import time
import itertools
import collections
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shutil import copyfile

# Make sure we take PBC into account
mda.core.periodic = True


def dim2lattice(x, y, z, alpha=90, beta=90, gamma=90):
    """Convert dimensions (lengths/angles) to lattice matrix"""
    cosa = np.cos( np.pi * alpha / 180 )
    cosb = np.cos( np.pi * beta / 180 )
    cosg = np.cos( np.pi * gamma / 180 )
    sing = np.sin( np.pi * gamma / 180 )

    zx = z * cosb
    zy = z * ( cosa - cosb * cosg ) / sing
    zz = np.sqrt( z**2 - zx**2 - zy**2 )

    return np.array([x, 0, 0, y * cosg, y * sing, 0, zx, zy, zz]).reshape((3,3))


def positive_linear_blur(array, box):
    """
    Perform simple line blurring of a matrix by rolling
    one time in positive x, y, and z directions, correcting
    for non-rectangular PBC.
    """
    blurred = np.copy(array)
    shift = 1
    # The box matrix is triangular
    
    # ... so a roll over x is just okay
    blurred += np.roll(blurred, shift, axis=0)
    
    # ... but a roll over y may have an x-shift
    #
    xshift = shift * box[1, 0]
    rolled = np.roll(array,  shift, 1)
    rolled[:, 0, :] = np.roll(rolled[:, 0, :], -xshift, 0)
    blurred += rolled
    
    # .. and a roll over z may have an x- and a y-shift
    #
    xyshift = shift * box[2, :2]
    rolled = np.roll(array,  shift, 2)
    rolled[:, :, 0] = np.roll(rolled[:, :, 0], xyshift, (0, 1))
    blurred += rolled

    return blurred


def linear_blur(array, box, span, inplace=True):
    """
    Perform linear blurring of an array by rolling
    over x, y, and z directions, for each value 
    up to span. If inplace is True, the rolled 
    array is always the target array, causing
    a full blur.
    """

    blurred = np.copy(array)

    if inplace:
        other = blurred
    else:
        other = array
        
    for shift in range(1, span+1):
        # The box matrix is triangular
        
        # ... so a roll over x is just okay
        blurred += np.roll(other, shift, axis=0)
        blurred += np.roll(other, -shift, axis=0)
        
        # ... but a roll over y may have an x-shift
        #
        xshift = shift * box[1, 0]
        rolled = np.roll(other,  shift, 1)
        rolled[:, 0, :] = np.roll(rolled[:, 0, :], -xshift, 0)
        blurred += rolled
        #
        rolled = np.roll(other, -shift, 1)
        rolled[:, -1, :] = np.roll(rolled[:, -1, :], xshift, 0)
        blurred += rolled
        
        # .. and a roll over z may have an x- and a y-shift
        #
        xyshift = shift * box[2, :2]
        rolled = np.roll(other,  shift, 2)
        rolled[:, :, 0] = np.roll(rolled[:, :, 0], xyshift, (0, 1))
        blurred += rolled
        #
        rolled = np.roll(other, -shift, 2)
        rolled[:, :, -1] = np.roll(rolled[:, :, -1], -xyshift, (0, 1))
        blurred += rolled
    return blurred

            
def blur(array, box, span):
    if span == -1:
        blurred = positive_linear_blur(array, box)
    elif span == 0:
        blurred = linear_blur(array, box, -1, inplace=False)
    else:
        blurred = linear_blur(array, box, span)
        #TODO This does not exist.
        #blurred = full_linear_blur(array, box, span)
    return blurred


def contour(array, box, span, inv=True):
    if inv:
        array = ~array
    return blur(array, box, span).astype(bool) ^ array
        

def voxelate_atomgroup(atomgroup, resolution, hyperres=False, max_offset=0.05):
    box = dim2lattice(*atomgroup.dimensions)
    # The 10 is for going from nm to Angstrom
    nbox = (box / (10 * resolution)).round().astype(int) # boxels
    unit = np.linalg.inv(nbox) @ box                     # voxel shape
    error = unit - 10 * resolution * np.eye(3)           # error: deviation from cubic
    deviation = (0.1 * (unit**2).sum(axis=1)**0.5 - resolution) / resolution

    if (np.abs(deviation) > max_offset).any():
        raise ValueError(
            'A scaling artifact has occured of more than {}% '
            'deviation from the target resolution in frame {} was '
            'detected. You could consider increasing the '
            'resolution.'.format(max_offset,
            atomgroup.universe.trajectory.frame)
        )

    transform = np.linalg.inv(box) @ nbox                 # transformation to voxel indices
    voxels = atomgroup.positions @ transform         
    if hyperres:
        # Blur coordinates 
        neighbors = hyperres * (np.mgrid[-1:2, -1:2, -1:2]).T.reshape((1, -1, 3))
        voxels = (voxels[:, None, :] + neighbors).reshape((-1, 3))
    voxels = voxels.astype(int)
        
    # Put everything in brick at origin
    for dim in (2, 1, 0):
        shifts = voxels[:, dim] // nbox[dim, dim]
        voxels -= shifts[:, None] * nbox[dim, :]
        
    return voxels, nbox


def gen_explicit_matrix(atomgroup, resolution=1, hyperres=False, max_offset=0.05):
    """
    Takes an atomgroup and bins it as close to the resolution as
    possible. If the offset of the actual resolution in at least one
    dimension is more than by default 5%, the function will stop and
    return an error specifying the actual offset in all dimensions
    plus the frame in which the mapping error occured.
    
    Returns
    (array) 3d boolean with True for occupied bins
    (dictionary) atom2voxel mapping

    """

    voxels, nbox = voxelate_atomgroup(atomgroup, resolution, hyperres, max_offset=max_offset)
    # Using np.unique gives a small performance hit.
    # Might still be necessary with large coordinate sets?
    # unique = np.unique(voxels, axis=0)
    x, y, z = voxels.T
    explicit = np.zeros(np.diagonal(nbox), dtype=bool)
    explicit[x, y, z] = True

    # generating the mapping dictionary
    voxel2atom = collections.defaultdict(list)
    # atom index starts from 0 here and is the index in the array, not the
    #  selection atom index in atom_select (these start from 1)
    if hyperres:
        indices = np.repeat(atomgroup.ix, 27)
    else:
        indices = atomgroup.ix
    for idx, voxel in zip(indices, voxels):
        voxel2atom[tuple(voxel)].append(idx)
        
    return explicit, voxel2atom, nbox


# def gen_explicit_matrix_multiframe(atomgroup, resolution=1,
#                                    max_offset=0.05, frames=0, hyper_res=False):
#     """
#     Tries to add multiple explicit matrices and mappings together. This might
#     be usefull for clustering at higher matrix resolution (lower than 1 nm).
#     Frames is used to smear over extra consecutive frames (0 is no smearing
#     over time). Hyper_res is used to smear the data points over half the
#     resolution. Hyper_res and frame smearing is exclusive for now.
    
#     Returns
#     (array) 3d boolean with True for occupied bins
#     (dictionary) atom2voxel mapping
#     """

#     if hyper_res:
#         return gen_explicit_matrix(
#             atomgroup, resolution, hyper_res, max_offset
#         )

#     # Starting the current frame
#     current_frame = atomgroup.universe.trajectory.frame
#     explicit_matrix, voxel2atom, nbox = gen_explicit_matrix(atomgroup, resolution, 
#                                                             False, max_offset)
    
#     # Try to stack the densities, but could fail due to voxel amount mismatch
#     #  due to pressure coupling and box deformations.
#     # Preventing index errors.
#     max_frame = len(atomgroup.universe.trajectory) - 1
#     end_frame = current_frame + frames
    
#     # actual expansion for multiframe smearing
#     frame = current_frame + 1
#     while frame <= max_frame and frame < end_frame:
#         atomgroup.universe.trajectory[frame]
#         # Smearing the positions for hyper res.

#         temp_explicit_matrix, temp_voxel2atom, nbox = gen_explicit_matrix(
#             atomgroup, resolution, max_offset
#         )
#         try:
#             explicit_matrix += temp_explicit_matrix
#             #voxel2atom = {**voxel2atom, **temp_voxel2atom}
#         except ValueError:
#             #TODO testing
#             #print('There was a mismerge.')
#             pass
#         frame += 1
#     # Set the active frame back to the current frame    
#     atomgroup.universe.trajectory[current_frame]
    
#     return explicit_matrix, voxel2atom, nbox


def voxels2atomgroup(voxels, voxel2atom, atomgroup):
    """
    Converts the voxels in a voxel list back to an atomgroup.
    
    Takes a voxel list and uses the voxel2atom mapping with respect to the
    atomgroup.universe to generate a corresponding atomgroup with the voxel 
    list. This is the inverse of gen_explicit_matrix.
    
    Returns an atomgroup.
    """
    # It is not important that every index only occurs onec,
    # as long as each atom is only selected once.
    indices = { idx for v in voxels for idx in voxel2atom[tuple(v)] }
    return atomgroup.universe.atoms[list(indices)]


def clusters2atomgroups(clusters, voxel2atom, atomgroup):
    """
    Converts the cluster in voxel space to an atomgroup.
    
    Clusters is a dictionary with the cluster ids as keys and the voxel_lists
    as values.
    
    Returns a list of atomgroups.
    """
    return [
        voxels2atomgroup(v, voxel2atom, atomgroup)
        for c, v in clusters.items()
    ]


def find_neighbours(position, box, span=1):
    """
    Uses the position to generate a a box width size span around the position.
    Taking cubic PBC into account.
    """
    s = slice(-span, span+1)
    neighbours = np.mgrid[s, s, s].T.reshape((-1, 3)) + position
    
    if all(m > span-1 for m in position) and all(position < np.diagonal(box) - span):
        # ... then we are done already
        return [ tuple(v) for v in neighbours ]

    # Put everything in brick at origin
    for dim in (2, 1, 0):
        shifts = neighbours[:, dim] // box[dim, dim]
        neighbours -= shifts[:, None] * box[dim, :]

    return [ tuple(v) for v in neighbours ]


#@profile
def non_clustered_atomgroup(atomgroup, cluster_array):
    """
    Takes an atomgroup and a cluster array of same size (2D array) and returns 
    an atomgroup containing all atoms in the atomgroup which have cluster 0 
    assigned in the cluster array.
    """
    # Find the non-clustered indices only for the headgroups with respect to 
    #  the headgroup indices
    non_clustered_indices = cluster_array[atomgroup.ix] == 0
    
    # Obtaining the atomgroup for which we have to perform a neighbourhood 
    #  search with a certain cutoff. This returns an atomgroup.
    non_clustered_atomgroup = atomgroup[non_clustered_indices] 
    
    return non_clustered_atomgroup


#@profile
def find_key_with_max_value(dictionary):
    """
    Takes a dictionary and returns the key with the highest value. 
    If there is no unique highest value, 0 is returned.
    
    EXPERIMENTAL WITH THE AT LEAST 5 BIGGER MAX VALUE.
    """
    #TODO PART OF TESTING THE MINIMUM REQUIREMENT
    min_diff = 0
    # convert the dictionary items to list so they are fixed
    values = list(dictionary.values())
    keys = list(dictionary.keys())
    # Find the maximum value
    max_value = dictionary[keys[values.index(max(values))]]
    # Find all keys with the maximum value
    max_keys = [keys for keys, values in dictionary.items() if 
                values == max_value]
    # Only return a hit if the maximum value is unique
    if len(max_keys) == 1:
        #TODO Maybe implement this in a nice manner I have to think about this
        for value in values:
            if max_value == value:
                continue
            if max_value <= value + min_diff:
                break
        else:
            return max_keys[0]
    
    return 0

#@profile    
def find_dominant_neighbour_cluster(ref_atomgroup, query_atomgroup, cutoff, 
                                    cluster_array, possible_clusters):
    """
    Uses the query_atomgroup to perform a neighbour search and find the 
    dominant cluster (if there is one) for a the residues of the atoms in 
    the query_atomgroup. If there is no unique dominant cluster, the function 
    returns 0. If there are cluster_array indices to alter, it will return 
    the query_atomgroup. The possible cluster can be given to prevent 
    searching for it in each iteration.
    """ 
    # Setting the reference atomgroup for the search, excluding the 
    #  self particles.
    ref = mda.lib.NeighborSearch.AtomNeighborSearch(
        ref_atomgroup - query_atomgroup, ref_atomgroup.dimensions,
        )
    # Performing the neighbourhood search with a cutoff of 10 Angstrom
    hits = ref.search(query_atomgroup, cutoff, 'A') # A is for Angstrom
    # Obtaining the cluster value for each hit if there is a hit at all.
    if len(hits) > 0:
        hit_clusters = cluster_array[hits.ix]
    else: 
        return 0
    
    # Counting the cluster prevalence in the hit_clusters
    cluster_count = {}
    for possible_cluster in possible_clusters:
        cluster_count[possible_cluster] = 0
    for hit_cluster in hit_clusters:
        cluster_count[hit_cluster] += 1
        
    # Obtaining the dominant cluster (if there is one)
    dominant_cluster = find_key_with_max_value(cluster_count)
    
    # Mapping the query headgroup complete residue indices with respect to the 
    #  cluster array and their dominant cluster. 
    if dominant_cluster != 0:
        return [query_atomgroup, dominant_cluster]
    # Return 0 is nothing was changed
    return 0

#@profile
def force_clustering(ref_atomgroup, cutoff, cluster_array, possible_clusters):
    """
    Assigns a clusters to all non clustered residues with a shared atom in the 
    query_atomgroup in the cluster array. The cluster ID is only changed if 
    there is a dominant cluster around the atom in the query_atomgroup within 
    a certain cutoff (no unique max --> no reassignment). The changes to the 
    cluster array are made in place, the funtion will return an empty list if
    it performed no alterations in the cluster_array and a list containing all
    clustered atomgroups, it also returns the leftover atoms. 
    
    All changes are made at once to make the fairest dominant cluster 
    assignment and prevent order dependency for the dominant cluster. An array
    of possble cluster can be specified to prevent recalculation
    by putting np.unique(cluster_array) at the position of possible clusters.
    
    #TODO I should also make the query_atomgroup an input for it can take
    a long time to compute for large systems and since we keep track of our
    changes, we should be able to update it using a change log apporach.
    """
    query_atomgroup = non_clustered_atomgroup(ref_atomgroup, cluster_array)
    query_residuegroup = query_atomgroup.residues
    
    # Try to make the changes and either still return an empy list or a list 
    #  of (atomgroup, dominant_cluster) or simply an empty list for a 
    #  failed case.
    changes = []
    leftovers = 0
    for residue in query_residuegroup:
        active_atoms = residue.atoms & ref_atomgroup
        temp_changes = find_dominant_neighbour_cluster(
            ref_atomgroup, active_atoms, cutoff, 
            cluster_array, possible_clusters,
            )
        # Only accept the change, if it returned non zero (the ouput for no 
        #  change). Else add them to leftovers
        if temp_changes == 0:
            leftovers += 1
        else:
            temp_changes[0] = temp_changes[0].residues.atoms
            changes.append(temp_changes)
    # Altering the cluster assignment in the cluster array for non 0 changes
    for change in changes:
            cluster_array[change[0].ix] = change[1]
    return changes, leftovers

#@profile
def iterative_force_clustering(ref_atomgroup, cutoff, cluster_array, 
                               possible_clusters, max_cutoff=20, max_stop=1, 
                               cutoff_increment=5, verbose=False):
    """
    This performs an iterative force_clustering and stops when there is 
    nothing changed. It returns a tuple(3) of occupation of each cluster per 
    iteration in a dict if verbose is set to true, as well as the dynamic 
    cutoffs list and the leftovers list. If verbose is not set, it will return 
    the iteration depth.
    """
    if verbose:
        unique, counts = np.unique(cluster_array, return_counts=True)
        output_dict = {}
        for idx, cluster in enumerate(unique):
            output_dict[cluster] = [counts[idx]]
    
    # Setting some initial parameters changes in this sense reflects the 
    #  amouot of groups that where considered by the algortym, this is due to
    #  the fact that a non-changed group returns a 0 and this also ends up in 
    #  the list. If no group are left which return either 0 or an atomgroup, 
    #  the changes are 0.
    start_cutoff = cutoff
    cutoffs = []
    changes = []
    old_change = 0
    change = -1
    leftovers = -1
    counter = 0
    stop = 0
    
    while stop != max_stop and leftovers != 0:
        # Perform force clustering for each non clustered residue in the 
        #  ref_atomgroup.
        changed_cluster_array_indices, leftovers = force_clustering(
            ref_atomgroup, cutoff, 
            cluster_array, possible_clusters,
            )
        # Some bookkeeping for proper quality control
        if verbose:
            unique, counts = np.unique(cluster_array, return_counts=True)
            temp_dict = dict(zip(unique, counts))
            cutoffs.append(cutoff)
            for key in temp_dict:
                output_dict[key].append(temp_dict[key])
        
        # Calculate the amount of changes
        change = len(changed_cluster_array_indices)
        
        # If nothing has changed, increment the cutoff and start the death 
        #  counter 
        if change == old_change:
            if (cutoff + cutoff_increment) <= max_cutoff:
                cutoff += cutoff_increment
            stop += 1
        # If something has changed, move to intial cutoff and reset death 
        #  counter
        else:
            cutoff = start_cutoff
            stop = 0
            
        # Updating change
        old_change = change
        changes.append(change)
        counter += 1
    
    if verbose:
        return output_dict, cutoffs, changes
    else:
        return counter

# The 3d example of a very clear more set oriented neighbour clustering
#@profile
def set_clustering(explicit_matrix, box, exclusion_mask=False, span=1, 
                   verbose=False):
    """
    A set oriented neighbour voxel clustering.
    
    The explicit matrix should be a 3d boolean array. The exclusion mask has
    to have the same dimensions as the explicit matrix and should also be 
    a boolean array. The exclusion mask acts as a dead zone for clustering. 
    The span is used to expand clustering to the first N neighbours
    in voxel space. Verbose can be used for more information during the 
    clustering.
    
    Returns a dictionary of clusters with a list of voxels per cluster.    
    """
    # Obtaining a set for all occupied voxels
    # starting the timer for generating the set (verbose)
    if verbose:
        start = time.time()
    # removing the exclusion mask from the occupied voxel matrix
    if exclusion_mask is not False:
        explicit_matrix[exclusion_mask == True] = False
    # finding all hits (occupied voxels) and
    # adding each hit to the set as a tuple
    pset = { pos for pos in zip(*np.where(explicit_matrix)) }
    if verbose:
        print('There are {} points to cluster.'.format(len(positions_set)))
    # stop timer for making the hits set (verbose)
    stop = time.time()
    if verbose:
        print('It took {} to make the set.'.format(stop-start))

    # The clustering scales linear and millions of points can be achieved in 
    #  the minute range.
    # beginning the timer for clustering (verbose)
    if verbose:
        start = time.time()
    # the starting cluster
    current_cluster = 1
    # output dictionary containing a list of hits per cluster
    clusters = {}
    # as long as there are hits
    while pset:
        clusters[current_cluster] = []
        active = {pset.pop()}
        while active:
            around = { n for v in active for n in find_neighbours(v, box, span) }
            clusters[current_cluster].extend(active)
            active = { n for n in around if n in pset and not pset.remove(n) }
        current_cluster += 1
    # stopping the timer for clustering (verbose)
    if verbose:
        stop = time.time()

    if verbose:
        print('It took {} to cluster {} clusters \
with a total of {} points'.format(stop-start, len(clusters), len(positions)))
            
    return clusters


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


def vmd_visualization_single_frame(template_file, clusters):
    """
    A small script to generate a vmd visualization to check the clustering
    manually.
    """
    alphabet = 'abcdefghijklmnopqrstuvwxyz'.upper()
    
    target_path = template_file.split('/')
    if len(target_path) > 1:
        target_path = '/'.join(target_path[:-1])
    else:
        target_path = '.'
            
    target_file = template_file.split('/')[-1]
    
    target_file = target_file.split('.')
    if len(target_file) > 1:
        target_file = '.'.join(target_file[:-1])
    else: 
        target_file = str(target_file)
    target_file += '_clustered.vmd'
    
    full_target_path = target_path + '/' + target_file
    print(full_target_path)
    
    copyfile(template_file,
             full_target_path)

    print(np.unique(clusters))
    with open(full_target_path, 'a') as f:
        for cluster in np.unique(clusters):
            f.write('\n\nset cluster{} '
                    '[atomselect top "index '.format(cluster))
            atom_indices = np.asarray(np.where(clusters == cluster))[0]
            print(atom_indices)
            for idx, element in enumerate(atom_indices):
                if idx % 10 == 0:
                    f.write('\\\n')
                f.write('{} '.format(element))
            f.write('"]\n\n')
            f.write('$cluster{0} set chain {1}'.format(
                    cluster, 
                    alphabet[((cluster-1) % 26)],
                    ))
