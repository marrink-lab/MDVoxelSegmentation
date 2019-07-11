###############################################################################
########### MDVoxel Clustering Beta input file Leaflets #######################
###############################################################################
###   All lines indicated as mandatory contain a variable which must be set. 
### The other lines are just a another way of obtaining the needed selections
### witout typing too much. Although they are also there to give an example
### on how to tackle the more challenging selections you might want to make.
###
### For now it is hard coded that this file should be named:
### clustering_input.py
###############################################################################
###############################################################################


### Input files (mandatory).
# You should be able to create a complete universe with residue information.
tpr = 'md.tpr'
# The actual xtc or gro
xtc = 'centered_whole_trajout.xtc'


### Output file (mandatory).
output_file = 'clusters'


### Full lipids (mandatory).
# Specifying the lipid residue names and concatenating the list.
lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS',
          'LYPC', 'LYPS', 'OA-1', 'OAOH']
lipids_query = ' '.join(lipids)
# Specifying the mda selection.
lipids_selection_query = 'resname {}'.format(lipids_query) # (mandatory)


### Tail of lipids (mandatory).
# Specify the tail beads to take into account and concatenating the list.
tails = ['C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A',
         'D3B', 'D2A', 'D2B', 'C2A', 'C2B']
tails_query = ' '.join(tails)
# Specifying the mda selection.
tails_selection_query = 'name {}'.format(tails_query) # (mandatory)


### Headgroups selection (mandatory).
# Specify the headgroups used to mask the lipid tail density.
headgroups = ['PO4', 'NC3', 'CNO', 'NH3', 'TAP', 'GL1', 'GL2']
headgroups_query = ' '.join(headgroups)
# Specifying the mda selection, set to false if not used (mandatory).
headgroups_selection_query = 'name {}'.format(headgroups_query) # (mandatory)
#headgroups_selection_query = False # Turns off the headgroups selection.


### Exclusions selection (mandatory).
# Specifying the exclusion selection and concatenating the list.
#exclusions_selection = ['BB', 'SC1', 'SC2','SC3', 'SC4']
#exclusions_query = ' '.join(exclusions_selection)
# Specifying the mda selection, set to False to turn off exclusions. 
#exclusions_selection_query = 'name {}'.format(exclusions_query) # (mandatory)
# Turning off exclusions.
exclusions_selection_query = False


### Threads for running the clustering (set to computer core count).
# reduce if you run out of memory. The amound of threads should be
# smaller than the number of frames in question. (mandatory).
threads = 12


### Setting the aimed voxel size in nm (mandatory).
resolution = 0.5


### Hyper res blurs the positions for voxel mapping over
# half the resolution. This is useful for high res CG clustering.
# such as the tight packing of an inverted hexagonal phase (mandatory).
hyper_res = True


### Setting foce lipid clustering by headgroup neighbourinvg clusters (mandatory).
force = True
force_cutoff = 20
force_info = False


### Settting the scatterplot output (mandatory).
plotting = False


### Setting the reducting in points for plotting (mandatory).
reduce_points = 19


### Setting the amount of prints and test plots (mandatory).
verbose = False


### start/stop/skip set stop to None to go to end of file (mandatory).
start_frame = 0
stop_frame = None
skip = 1


### Extra frames for averaging over time (experimental, mandatory)
frames = 0


### Actually running the leaflet clustering my executing this file with python3.
if __name__=='__main__':
    import mdvoxelclustering as mdv
    # Run the position based leaflet clustering (outputs clusters.npy).	
    #mdv.leaflets.main()
    mdv.leaflets.main_threaded()
    # Run the clluster identity over time (outputs clusters_ordered.npy).
    mdv.settests.main()
