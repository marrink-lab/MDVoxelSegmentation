### Input files
# you should be able to create a complete universe (mandatory)
tpr = '../PROD/md.tpr'
# the actual xtc or gro (mandatory)
xtc = '../PROD/md.xtc'

### Output file
output_file = 'clusters'


### Tail of lipids
# specify the tail beads to take into account and concatenating the list
tails = ['D6A', 'D6B', 'C6A', 'C6B', 'C5A', 'C5B', 'D5A', 'D5B','C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A', 'D3B',]# 'D2A', 'D2B', 'C2A', 'C2B']
tails_query = ' '.join(tails)
# specifying the mda selection, this is mandatory to have!
tails_selection_query = 'name {} or (resname CHOL XHOL and name C1)'.format(tails_query)


### Headgroups selection
# specify the headgroups used to mask the lipid tail density
headgroups = ['COO', 'COOH', 'PO4', 'NC3', 'CNO', 'NH3', 'TAP', 'GL1', 'GL2', 'AM1', 'AM2', 'GM1', 'GM2', 'GM3', 'GM4', 'GM5', 'GM6', 'GM7', 'GM8', 'GM9', 'GM10', 'GM11', 'GM12, GM13', 'GM14', 'GM15', 'GM16', 'GM17']#, 'C1A', 'C1B', 'C2A', 'C2B', 'D1A', 'D1B', 'D2A', 'D2B']
headgroups_query = ' '.join(headgroups)
# specifying the mda selection, set to false if not used (mandatory)
headgroups_selection_query = 'name {} or (resname CHOL and name ROH) or (resname PAPI PIPI POP1 POP2 POP3 POPI PUPI and name C1 C2 C3 P1 P2 P3)'.format(headgroups_query)
#headgroups_selection_query = False


### Exclusions selection
# specifying the exclusion selection and concatenating the list
exclusions_selection = ['BB', 'SC1', 'SC2','SC3', 'SC4']
exclusions_query = ' '.join(exclusions_selection)
# specifying the mda selection, set to False to turn off exclusions. 
# This is mandatory to have!
exclusions_selection_query = 'name {}'.format(exclusions_query)
# Turning off exclusions
#exclusions_selection_query = False

### verbose for testing
verbose = False

### start stop and skip
start_frame = 0
# Use non to set complete trajector
stop_frame = None
skip = 1

# Setting the minimum size in particles for the cluster
min_cluster_size = 50

### Setting the aimed voxel size in nm (mandatory) als has the frames
# used for smearing in time (0 is no time smearing).
resolution = 0.5
frames = 0
hyper_res = True

### force clustering of masked residues
force = True
force_cutoff = 20
force_info = False

### Threads
threads = 10

### Setting the trajectory skip size
skip = 1

### Settting the scatterplot output
plotting = False

### Setting the reducting in points for plotting
reduce_points = 1


if __name__=='__main__':
    import mdvoxelclustering as mdv
    mdv.leaflets.main_threaded()
    mdv.settests.main()
    mdv.antidisco.main()
    mdv.graphing.main()
