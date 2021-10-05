===============================
MDVoxelSegmentation
===============================
Using neighbor segmentation in voxelspace for fast and consistent spatial and temporal segmentation.

This software has been developed to allow for a higher selection syntax than atom and or residue index, such as abstract complex particles (e.g. lipid monolayers). MDVoxelSegmentation combines multiple layers of neighbor segmentation with a voxel mask, and makes tracking segments over time possible at high quality.

This code was published in the following `article <https://doi.org/10.1021/acs.jctc.1c00446>`_.

* Open software: Apache 2 license

.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png

Features
--------
* v1.0 is a stable build and segmentation should be of high quality. The code is usable for high throughput with minimal effort with optimization. By default no more than a GRO and XTC (or equivalent) are required for successful segmentation of Martini systems.
* Voxel based neighbour segmentation under all periodic boundary conditions
* Fast contour segmentation
* Compatible with most MD file formats due to its tight link to MDAnalysis
* Consistent segmentation over time on trajectories
* Compatible with VMD using standard visualization files (python compiled VMD)
* A wide range of examples (bottom of this readme)
* Membrane leaflet assignment of lipids of most topologies
    - Bilayers
    - Vesicles
    - (Inverted) hexagonal
    - Cubic
    - Membrane tethers
    - Complex lipids formulations including cholesterol
    - Proteins
    - Up to millions of beads per frame (possibly billions)
    
Instructions
--------
Installation
************
Fork the development branch using:

:code:`git clone https://github.com/BartBruininks/MDVoxelSegmentation.git`

Move into the cloned folder and type:

:code:`pip install -e .`

Add an alias in your :code:`~/.bashrc` and restart/resource, finally run mdvseg with the help flag:

:code:`alias mdvseg="python /path/to/mdvsegmentation_package/mdvoxelsegmentation/do_segmentation.py"` (add to :code:`~/.bashrc` and make sure the path is pointing to the correct folder)

:code:`source ~/.bashrc` (run in terminal)

:code:`mdvseg -h` (run in terminal)

Basic Segmentation
***************
To perform default segmentation on a GRO and XTC file containing a coarse grain Martini system, you have to specify the GRO and XTC file. The final segmentation assignment will be written to :code:`clusters_ordered.npy`. This file can be used using numpy in python to perform the required analysis. 

:code:`mdvseg -f path_to_your.gro -x path_to_your.xtc`

MDAnalysis will probably throw some warnings stating that it cannot estimate the masses for you coarse grain particles. This cannot be suppressed but it harmless. A useful graph of your segmentation can be made by running in the output folder:

:code:`python /path/to/mdvoxelsegmentation_package/mdvoxelsegmentation/plotting.py`

The graph is created by default and makes use of matplotlib. Therefore adjusting the plotting script to make exactly what you need should be rather straight forward.

It is mainly the force segmentation flag (-fs) and it associated recursion depth (-rd) which have a big impact on performance. Turning force segmentation off (-lg 0) is fine if perfect final quality is not needed. By default `mdvseg` generates its own `selections.inp` which should cover all basic lipids in Martini. However, if some definitions are missing, you can always manually add them to the `selections.inp` (mdvseg does NOT overwrite an already present `selection.inp`).

VMD visualization
******************
For visualization with VMD you need to make sure that you have a version of VMD compiled against a python version supporting numpy. A compatible VMD compilation will be distributed in the future, for now you have to figure this out yourself, I did ask the developers to support anyone asking for such compilation and they said yes! So just send an e-mail to the VMD mailing list if you would need it. If you have the right flavor of VMD, all you need to type next is:

:code:`vmd -e vmd_clusters_visualization.vmd`

The first 32 segments will automatically be assigned a color and material/style. They can be used to make selections using 'user your_segment' in the VMD selection syntax. These representations should automatically be set to update every frame. By typing 'hide' in the VMD terminal, we can easily turn off all segment representations. Segment 0 always contains everything which was not assigned a segment and is hidden. 'user' 32 always shows segment 32 to 1000, to show all segments which might have a very high index. The downside is that all segments from 32 onward have the same color.

Useful things to know
*********************
Using MDVoxelSegmentation on coarse grain Martini lipid/protein systems should work without needing much prior knowledge. However, to make the most out of the created :code:`clusters_ordered.npy` it is useful to know some python (numpy, MDAnalysis, Matplotlib). If you are working with atomistic systems and have to specify your own headgroups/linkers/tails, you need to known what the relevant names are from your PDB/GRO and make your own selection entries in the :code:`selections.inp`. The :code:`selections.inp` uses the MDAnalysis selection syntax (very close to the VMD selection syntax). Below are some basic lines of code to help you on your way with using the segmentation data. First we will give an example for some basic plotting, followed by an example for an atomistic CHARMM :code:`selections.inp` for DOPE lipids.

*A basic python example to plot the number of segments over time*

.. code-block:: python

    ## Importing numpy and matplotlib.
    import numpy as np
    import matplotlib.pyplot as plt

    ## Loading the segmentation data.
    segments_over_time = np.load('clusters_ordered.npy')

    ## Calculating the amount of segments in each frame.
    # Make an empty array which has one int32 for every frame.
    segments_per_frame = np.zeros(segments_over_time.shape[0], dtype='int32') 
    # Fill the array with the amount of non-zero segments in each frame.
    segments_per_frame[:] = [len(np.unique(frame)) - 1 for frame in segments_over_time]

    ## Plotting the results.
    # Making an empty plot.
    fig, ax = plt.subplots()
    # Adding required data to plot.
    ax.plot(segments_per_frame)
    # Set ticks to a sensible regime.
    start, end = [round(limit) for limit in ax.get_ylim()]
    ax.yaxis.set_ticks(np.arange(start, end + 1, 1))
    # Add labels to axes.
    ax.set_xlabel('Frame count')
    ax.set_ylabel('Number of segments')
    # Save the plot.
    fig.savefig('amount_of_segments_over_time.png', dpi=300)
    # Usually people don't like it if you pop windows, however
    #  if you would like to automatically show the result uncomment
    #  the following line.
    #fig.show()

*An atomistic segmentation example for DOPE lipids with the CHARMM force field*

.. code-block:: python

    $ vi selections.inp
    ## Create an empty `selections.inp` and add the following lines, the selection 
    ##  syntax should always be one line and directly follow its header description.

    # It is not bad to include the linkers also in the headgroups, but this is often 
    #  not important.
    [charmm_heads]
    (name N P C12 C11 O11 O12 O13 O14)

    [charmm_linkers]
    (name C1 C2 O21 C21 C3 O31 C31)

    # Not adding the first carbons of the tail can improve performance, but again, 
    #  this is usually not important.
    [charmm_tails]
    (name C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C32 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C317 C318)


    ## Run the mdvseg, hyper resolution can be turned off for there are more than
    ##  enough particles due to the atomistic resolution!
    $ mdvseg -f your.gro -x your.xtc -hg charmm_heads -lg charmm_linkers -tg charmm_tails -hres 0

Post some feedback in our issues
*********************************
We would love to see how you used MDVoxelSegmentation in your projects and are always happy to see cool systems and screenshots. Just upload them to an issue of this repository. There is also place for comments on usability and bugs. 

How to contribute
******************
If you are interested in joining this project after its initial release just post an issue, or better yet, send an email to bartbruininks@gmail.com. We are currently still filled with ideas to be implemented and all hands are welcome. Some of the open topics are in the issues. Anybody who contributes for a fair share will off course be included in future publications.

Examples
---------
.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png
**Segmentation of the inverted hexagonal phase with four inner channels connected to a bilayer with a fusion stalk**

Inside the channels is a fragment of dsDNA. The leaflet segmentation was performed using a resolution of 0.5 and hyperesolution turned on. This to allow for the correct segmentation of the tight geometry of the channels in coarse grain data (Martini, we used hyper resolution for all CG data!), also force segmentation was turned on to have (almost?) every lipid assigned up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61180812-f9b20780-a61c-11e9-838f-f42e54133669.png
**Leaflet segmentation of a complex plasmamembrane tether**

The two leaflets of the plasmamembrane are clearly assigned correctly and depicted as a transparent surface. The cholesterol inside the two leaflets is drawn in VDW spheres and their headgroups have a slightly altering color. All cholesterol seems to be assigned correctly. Segmentation was performed with a 0.5 nm resolution and iterative forced segmentation to assign the diving cholesterol up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/75271704-e7c45400-57fc-11ea-896a-60f0e2718f0d.png
**Leaflet segmentation of a plasma membrane including multiple proteins**

Less than 30 lipids remain unassigned of the roughly 1 million present. The leaflet assignment seemed to have worked correctly. For segmentation a resolution of 0.5 nm and iterative forced segmentation within 2 nm was used. The protein was used as exclusion to prevent them acting as pores in our segmentation. In total 1.3 millions beads were segmented in 30 minutes on a desktop. Mainly the force segmentation to assign all diving leaflets took a while. Keep in mind that this will never change the amount of segments present, so forced segmentation could be skipped in many situations. This is also only making use of a single core (a single frame cannot be hypterthreaded in the current code).

.. image:: https://user-images.githubusercontent.com/1488903/75272814-e3009f80-57fe-11ea-868d-29b1bd126c7a.png
**A collection of notoriously hard bilayer problems for segmentation**

For the cholesterol flip-flopping we use non iterative forced segmentation with a cutoff of 1.5 nm to act as a deadzone of 1 nm (A, B, C; recursion depth set to 1). We see that intercalating close contact leaflets do not cause faulty segmentation (D, E). Pores are also handled correctly and the minimum pore size at a resolution of 0.5 nm is actually the pore itself (F/G). If the pore is only a water channel, but the lipids do not reorient, its not considered a pore. Since the leaflets are not even continuous. In short we are able to detect all *toroidal* pores in a membrane. Water pores are a different game which we might solve in the future with a similar set based approach (ohh yhea we got something nice brewing, if only we had time :D).

.. image:: https://user-images.githubusercontent.com/1488903/75491447-4a148480-59b6-11ea-92ef-6faf0c646333.png
**Single frame toroidal and/or water pore detection in a bilayer**

A small glimps of what we are workin on with the pores. As you can see we can identify both toroidal (left) and water only pores (right). The frames were handpicked for we looked specifically for a toroidal and water pore. The expected end goal would be the consistent identification of all pores in membranes. Just as we do for leaflets. The pore tracking should be combinable with the leaflet identification, allowing for segmentation using the pores as exclusion mask. At the same time the pores would just have their own segmentation array which can be used for later analysis and visualization. This allows for leaflet identification, even in the presence of water and/or toroidal pores.

Credits
---------
Bart M. H. Bruininks, Albert Thie, Paulo C. T. de Souza, Tsjerk A. Wassenaar, Shirin Faraji & Siewert J. Marrink

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
