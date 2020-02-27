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

This software has been developed to allow for a higher selection syntax than atom and or residue index, such as abstract complex particles (e.g. lipid monolayers). MDVoxelClustering combines nieghbour clustering with a voxel mask of set resolution (default 0.5 nm). Forced-segmentation is turned on by default and assigns particles to clusers if they were unassigned. This happens within a 2 nm cutoff in an iterative manner. Finally there is a minimum cluster size of 50 particles to prevent clutter. For now the software has mainly been tested on CG Martini lipid assemblies with or without embedded proteins, though it probably works even better on atomistic simulations due to the higher amount of atoms that can be binned. If working with atomistic systems, you can probably turn off hyper-resolution in the input file.

.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png

* Open software: Apache 2 license

Features
--------
* v0.9 is a beta build and segmentation should be checked by eye. High throughput is not adviced yet. Exmple input files for CG lipid systems are supplied. They can be easily adapted to your system by quickly following the steps in the instructions.
* Voxel based neighbour clustering under cubic periodic boudaries
* Fast contour clustering
* Compatible with most MD file formats due to its tight link to MDAnalysis
* Consistent clustering over time on trajectories of any size
* Compatible with VMD using standard visualization files
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
Some short instructions on using the example files on a CG lipid containing system.

Installation
************
Fork the development branch using:

:code:`git clone [branch_clone_url]`

Then move into the cloned folder and type:

:code:`pip install -e .`

Setting up our input file for CG leaflet segmentation
******************************************************
Example files for leaflet clustering and VMD visualization can be found in the `example_inputs` folder.

The basic is that we set the correct path to 'our.tpr and 'our.gro'/'our.xtc' in the 'clustering_input.py' file. This is best done by making a copy of this file into an anlysis direcotry. The paths are set directly at the top of the file. We also need to set the amount of threads we want to use, by default this is set to 12. The trajectory should contain at least as many frames as the amount of threads we set. Multithreading for now is just chunking your trajectory in pieces and handling them independently. Therefore memory consumption roughly scales linear with the amount of threads.

Running the segmentation and creating sensible output
******************************************************
If 'mdvoxelclustering' is installed in our local python3 environment, we can run the clustering by executing the 'clustering_input.py' in our terminal:

:code:`python clustering_input.py`

In the bottom of the 'clustering_input.py' we can see why this works. The input_file itself runs the necessary python functions to perform leaflet segmentation and then improves the raw output by running a consistent id algorithm. Finally our segments will be plotted in 'clusters_over_time.png' as a graph and stored in 'clusters_ordered.npy' as a compressed numpy array. Some other files are created which are mainly used for plotting right now. You do not need to worry about those for now and just take a look at the created graph. Maybe you can already spot what is going on. The graph is a walking average of the cluster atom count, where grey lines indicate interactions between segments. The hue of the vertical line indicates the relative amount of particles in the system which participate in that interaction.

VMD visualization
******************
For futher visualization in VMD we need to add 'our.gro' and 'our.xtc' path to 'vmd_clusters_visualization.py'. These entries are somehwere halfway the file... (sorry). We also need to make sure that we have a version of VMD compiled against a python version supporting numpy. A compatible VMD compilation will be distributed in the future, for now you have to figure this out yourself, though I asked the developers to support anyone asking for such compilation and they said yes! So just send an e-mail to the VMD mailing list if you would need it. If we have the right falvour of VMD, all we need to type next is:

:code:`vmd -e vmd_clusters_visualization.vmd`

Our first 32 segments will automatically be assigned a color and material/style. They can be used to make selections using 'user your_segment' in the VMD selection syntax. These representations should automatically be set to update every frame. By typing 'hide' in the VMD terminal, we can easily turn off all segment representations. Segment 0 always contains everything which was not assigned a segment and is hidden. 'user' 32 always shows segment 32 to 1000, to show all segments which might have a very high index. The downside is that all segments from 32 onwards have the same color.

Post some feedback in our issues
*********************************
We would love to see how you used MDVoxelClustering in your projects and are always happy to see cool systems and screenshots. Just upload them to an issue of this repository. There is also place for comments on usability and bugs. 

How to contribute
******************
If you are interested in joining this project after its initial release just post an issue, or better yet, send an email to b.m.h.bruininks@gmail.com. We are currently still filled with ideas to be implemented and all hands are welcome. Some of the open topics are in the issues. Anybody who contributes for a fair share will off course be included in future publications.

Examples
---------
.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png
**Clustering of the inverted hexagonal phase with four inner channels connected to a bilayer with a fusion stalk**

Inside the channels is a fragment of dsDNA. The leaflet clustering was performed using a resolution of 0.5 and hyperesolution turned on. This to allow for the correct clustering of the tight geometry of the channels in coarse grain data (Martini, we used hyper resolution for all CG data!), also force clustering was turned on to have (almost?) every lipid assigned up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61180812-f9b20780-a61c-11e9-838f-f42e54133669.png
**Leaflet clustering of a complex plasmamembrane thether**

The two leaflets of the plasmamembrane are clearly assigned correctly and depicted as a transparent surface. The cholesterol inside the two leaflets is drawn in VDW spheres and their headgroups have a slightly altering colour. All cholesterol seems to be assigned correctly. Clustering was performed with a 0.5 nm resolution and iterative forced clustering to assign the diving cholesterol up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/75271704-e7c45400-57fc-11ea-896a-60f0e2718f0d.png
**Leaflet clustering of a plasma membrane including multiple proteins**

Less than 30 lipids remain unassigned of the roughly 1 million present. The leaflet assignment seemed to have worked correctly. For clustering a resolution of 0.5 nm and iterative forced clustering within 2 nm was used. The protein was used as exclusion to prevent them acting as pores in our segmentation. In total 1.3 millions beads were clustered in 30 minutes on a desktop. Mainly the force clustering to assign all diving leaflets took a while. Keep in mind that this well never change the amount of segments present, so forced clustering could be skipped in many situations. This is also only making use of a single core (a single frame cannot be hypterthreaded in the current code).

.. image:: https://user-images.githubusercontent.com/1488903/75272814-e3009f80-57fe-11ea-868d-29b1bd126c7a.png
**A collection of notoriously hard bilayer bilayer problems for segmentation**

For the cholesterol fli-flopping we use non iterative forced clustering (currently hard coded) with a cutoff of 1.5 nm to act as a deadzone of 1 nm (A, B, C). We see that intercalating close contact leaflets do not cause faulty segmentation (D, E). Pores are also handled correctly and the minimum pore size at a resolution of 0.5 nm is actually the pore itself (F/G). If the pore is only a water channel, but the lipids do not reorient, its not considered a pore. Since the leaflets are not even continous. In short we are able to detect all *toroidal* pores in a membrane. Water pores are a different game which we might solve in the future with a similar set based approach (ohh yhea we got something nice brewing, if only we had time :D).

.. image:: https://user-images.githubusercontent.com/1488903/75491447-4a148480-59b6-11ea-92ef-6faf0c646333.png
**Single frame toroidal and/or water pore detection in a bilayer**

A small glimps of what we are workin on with the pores. As you can see we can identify both toroidal (left) and water only pores (right). The frames were handpicked for we looked for specifically for a toroidal and water pore. The expected end goal would be the consistent identification of all pores in membranes. Just as we do for leaflets. The pore tracking should be combinable with the leaflet identification, allowing for segmentation using the pores as exclusion mask. At the same time the pores would just have their own segmentation array which can be used for later analysis and visualization. This allos for leaflet identification, even in the presence of water and/or toroidal pores.

Credits
---------
Bart M. H. Bruininks & Albert Thie

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
