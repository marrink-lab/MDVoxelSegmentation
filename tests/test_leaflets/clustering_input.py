###############################################################################
#################### MDVoxel Clustering Beta input file #######################
###############################################################################

###   All lines indicated as mandatory contain a variable which must be set. 
### The other lines are just a another way of obtaining the needed selections
### witout typing too much. Although they are also there to give an example
### on how to tackle the more challenging selections you might want to make.

### Input files (mandatory)
# you should be able to create a complete universe.
tpr = 'test_bilayer.tpr'
# the actual xtc or gro
xtc = 'test_bilayer.tpr'


### Output file (mandatory)
output_file = 'clusters'


#### Full lipids
## specifying the lipid residue names and concatenating the list
#lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS', 
#          'LYPC', 'LYPS']
#lipids_query = ' '.join(lipids)
## specifying the mda selection
#lipids_selection_query = 'resname {}'.format(lipids_query) # (mandatory)


### Tail of lipids
# specify the tail beads to take into account and concatenating the list
tails = ['C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A', 
         'D3B', 'D2A', 'D2B', 'C2A', 'C2B']
tails_query = ' '.join(tails)
# specifying the mda selection, this is mandatory to have!
tails_selection_query = 'name {}'.format(tails_query) # (mandatory)

### Headgroups selection
# specify the headgroups used to mask the lipid tail density
headgroups = ['PO4', 'NC3', 'CNO', 'NH3', 'GL1', 'GL2']
headgroups_query = ' '.join(headgroups)
# specifying the mda selection, set to false if not used (mandatory)
headgroups_selection_query = 'name {}'.format(headgroups_query) # (mandatory)
#headgroups_selection_query = False # Turns off the headgroups selection


### Exclusions selection
# specifying the exclusion selection and concatenating the list
exclusions_selection = ['BB', 'SC1', 'SC2','SC3', 'SC4']
exclusions_query = ' '.join(exclusions_selection)
# specifying the mda selection, set to False to turn off exclusions. 
exclusions_selection_query = 'name {}'.format(exclusions_query) # (mandatory)
# Turning off exclusions
exclusions_selection_query = False


### Setting the aimed voxel size in nm (mandatory)
resolution = 1


### Setting if you want to use the outer two contours for leaflets (mandatory)
expand = True # Does this even do anything anymore?


### Settting the scatterplot output (mandatory)
plotting = True


### Setting the reducting in points for plotting (mandatory)
reduce_points = 1


### Setting the amount of prints and test plots (mandatory)
verbose = False


### start/stop/skip (mandatory)
start_frame = 0
stop_frame = None
skip = 1
