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

This software has been developed to allow for a higher selection syntax than atom and or residue index. MDVoxelClustering combines nieghbour clustering with a voxel mask of set resolution (default 0.5 nm). Forced-segmentation is turned on by default and assigns particles to clusters if they were uncleastered within 2 nm in an iterative manner. Finally there is a minimum cluster size of 50 particles to prevent clutter. For now the software has mainly been tested on CG Martini lipid assemblies with or without embedded proteins, though it probably works even better on atomistic simulations due to the higher amount of atoms that can be binned. If working with atomistic systems, you can probably turn off hyper-resolution in the input file.

* Open software: Apache 2 license

Features
--------
* v0.9 is a beta build and visualization should be checked by eye. High throughput is not adviced yet. Exmple input files for CG lipid systems are supplied. They can be easily adapted to your system by quickly following the steps in the instructions.
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
    - Proteins
    - Up to millions of beads in seconds per frame
Instructions
--------
Some short instructions on using the example files on a CG lipid containing systems.

Installation
************
The '-e' is for editable install, which allows you to make changed in the library and python will actually check for those.

pip install -e .

Setting up our input file for CG leaflet segmentation
******************************************************
Example files for leaflet clustering and VMD visualization can be found in the `example_inputs` folder.

The basic is that we set the correct path to 'our.tpr and 'our.gro'/'our.xtc' in the 'clustering_input.py' file. This is best done by making a copy of this file into an anlysis direcotry. The paths are set directly at the top of the file. We also need to set the amount of threads we want to use, by default this is set to 12. The trajectory should contain at least as many frames as the amount of threads we set. Multithreading for now is just chunking your trajectory in pieces and handling them independently. Therefore memory consumption roughly scales linear with the amount of threads.

Running the segmentation and creating sensible output
******************************************************
If 'mdvoxelclustering' is installed in our local python3 environment, we can run the clustering by executing the 'clustering_input.py' in our terminal:

python clustering_input.py

In the bottom of the 'clustering_input.py' we can see why this works. The input_file itself runs the necessary python functions to perform leaflet segmentation and then improves the raw output by running a consistent id algorithm. Finally our segments will be plotted in 'clusters_over_time.png' as a graph and stored in 'clusters_ordered.npy' as a compressed numpy array. Some other files are created which are mainly used for plotting right now. You do not need to worry about those for now and just take a look at the created graph. Maybe you can already spot what is going on. The graph is a walking average of the cluster atom count, where grey lines indicate interactions between segments. The hue of the vertical line indicates the relative amount of particles in the system which participate in that interaction.

VMD visualization
******************
For futher visualization in VMD we need to add 'our.gro' and 'our.xtc' path to 'vmd_clusters_visualization.py'. These entries are somehwere halfway the file... (sorry). We also need to make sure that we have a version of VMD compiled against a python version supporting numpy. A compatible VMD compilation will be distributed in the future, for now you have to figure this our yourself, though I asked the developers to support anyone asking for such compilation and they said yes! So just send a mail to the VMD mailing list if you would need it. All we now need to type is:

vmd -e vmd_clusters_visualization.vmd

Our first 32 segments will automatically be assigned a color and material. They can be used to make selection using 'user your_segment' in the VMD selection syntax. By typing 'hide' in the VMD terminal, one can easily turn off all segment representations. Segment 0 always contains everything which was not assigned a segment. 'user' 32 always shows segment 32 to 1000, to show all segments which might have a very high index. The downside is that all segments from 32 onwards have the same color.

Post some feedback in our issues
*********************************
We would love to see how you used MDVoxelClustering in your projects and are always happy to see cool systems and screenshots. Just upload them to an issue of this repository. Of course there is also place for comments on usability and bugs. 

How to contribute
******************
If you are interested in joining this project after its initial release just post an issue, or better yet, send an email to b.m.h.bruininks@gmail.com. We are currently still filled with ideas to be implemented and all hands are welcome. Some of the open topics are in the issues. Anybody who contributes for a fair share will of course be included in future publications.

Examples
---------
.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png
**Clustering of the inverted hexagonal phase with four inner channels connected to a bilayer with a fusion stalk.**

Inside the channels is a fragment of dsDNA. The leaflet clustering was performed using a resolution of 0.5 and hyperesolution turned on. This to allow for the correct clustering of the tight geometry of the channels in coarse grain data (Martini), also force clustering was turned on to have (almost?) every lipid assigned up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61180812-f9b20780-a61c-11e9-838f-f42e54133669.png
**Leaflet clustering of a complex plasmamembrane thether.**

The two leaflets of the plasmamembrane are clearly assigned correctly and depicted as a transparent surface. The cholesterol inside the two leaflets is drawn in VDW spheres and their headgroups have a slightly altering colour. All cholesterol seems to be assigned correctly. Clustering was performed with a 1 nm resolution and forced clustering to assign (all?) the diving cholesterol up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61181667-b90cbb00-a629-11e9-9fc0-b2d52e4eaa93.png
**Leaflet clustering of a plasma membrane including multiple proteins.**

The issue described is not present anymore in the new clustering, however we are still working on the final figures for the paper. We can promise you this issue has been completely resolved and less than 30 lipids remain unassigned of the roughly 1 million present. The leaflet assignment seemed to have worked correctly, however, we do see some noise in places where we wouldn't expect it at first sight and this behaviour is to be further inspected. For clustering a resolution of 1 nm and forced clustering within 2 nm was used. The protein was used as exclusion to prevent flopping lipids next to proteins to intervene with the leaflet assignment. In total 1.3 millions beads were clustered in less than a minute on a desktop.

Credits
---------
Bart M. H. Bruininks & Albert Thie

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
