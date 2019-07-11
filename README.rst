===============================
MDVoxelClustering
===============================

.. image:: https://img.shields.io/pypi/v/mdvoxelclustering.svg
        :target: https://pypi.python.org/pypi/mdvoxelclustering

.. image:: https://img.shields.io/travis/BartBruininks/mdvoxelclustering.svg
        :target: https://travis-ci.org/BartBruininks/mdvoxelclustering

.. image:: https://readthedocs.org/projects/mdvoxelclustering/badge/?version=latest
        :target: https://readthedocs.org/projects/mdvoxelclustering/?badge=latest
        :alt: Documentation Status


Using neighbour clustering in voxelspace for fast and consistant spatial and temporal clustering.

This software has been developed to allow for a higher selection syntax than atom and or residue index. MDVoxelClustering combines nieghbour clustering with a voxel mask of set resolution (default 1 nm). For now the software has mainly been tested on CG Martini lipid assemblies with or without embedded proteins, though it probably works even better on atomistic simulations due to the higher amount of atoms that can be binned.

* Open software: Apache 2 license
* Documentation: https://mdvoxelclustering.readthedocs.org.

Features
--------
* This is an alpha build and should be used with extreme care!
* Voxel based neighbour clustering under cubic periodic boudaries
* Fast contour clustering
* Compatible with most MD file formats due to its tight link to MDAnalysis
* Consistent clustering over time on trajectories of any size
* Membrane leaflet recognition of lipids of most topologies
    - Bilayers
    - Vesicles
    - Inverted hexagonal phase
    - Membrane thethers
    - Complex lipids formulations including cholesterol
    - Up to millions of beads in seconds

```python
###############################################################################
########### MDVoxel Clustering Beta input file Leaflets #######################
###############################################################################

###   All lines indicated as mandatory contain a variable which must be set. 
### The other lines are just a another way of obtaining the needed selections
### witout typing too much. Although they are also there to give an example
### on how to tackle the more challenging selections you might want to make.

### Input files (mandatory)
# you should be able to create a complete universe.
tpr = 'your_reference.tpr'
# the actual xtc or gro
xtc = 'your_reference.xtc'


### Output file (mandatory)
output_file = 'clusters'


### Full lipids
# specifying the lipid residue names and concatenating the list
lipids = ['DLPC', 'DLPS', 'DOPE', 'DOTAP', 'DYPC', 'DYPS', 
          'LYPC', 'LYPS', 'OA-1', 'OAOH']
lipids_query = ' '.join(lipids)
# specifying the mda selection
lipids_selection_query = 'resname {}'.format(lipids_query) # (mandatory)


### Tail of lipids
# specify the tail beads to take into account and concatenating the list
tails = ['C4A', 'C4B', 'D4A', 'D4B', 'C3A', 'C3B', 'D3A', 
         'D3B', 'D2A', 'D2B', 'C2A', 'C2B']
tails_query = ' '.join(tails)
# specifying the mda selection, this is mandatory to have!
tails_selection_query = 'name {}'.format(tails_query) # (mandatory)

### Headgroups selection
# specify the headgroups used to mask the lipid tail density
headgroups = ['PO4', 'NC3', 'CNO', 'NH3', 'TAP', 'GL1', 'GL2']
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


### Threads for running the clustering (set to computer core count)
# reduce if you run out of memory. The amound of threads should be 
# smaller than the number of frames in question.
threads = 12


### Setting the aimed voxel size in nm (mandatory)
resolution = 0.5


### Hyper res blurs the positions for voxel mapping over
# half the resolution. This is useful for high res CG clustering.
hyper_res = True


### Setting foce lipid clustering by headgroup neighbourinvg clusters
force = True
force_cutoff = 20
force_info = False

### Settting the scatterplot output (mandatory)
plotting = False


### Setting the reducting in points for plotting (mandatory)
reduce_points = 19


### Setting the amount of prints and test plots (mandatory)
verbose = False


### start/stop/skip (mandatory) set stop to None to go to end of file.
start_frame = 0
stop_frame = None
skip = 1


### Extra frames for averaging
frames = 0

if __name__=='__main__':
    import mdvoxelclustering as mdv
    #mdv.leaflets.main()
    mdv.leaflets.main_threaded()
    mdv.settests.main()
```

Credits
---------

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
