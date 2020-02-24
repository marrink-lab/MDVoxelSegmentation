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
* v0.9 is a beta build and visualization should be checked by eye. High throughput is not adviced yet.
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

Examples
--------
Example files for leaflet clustering and VMD visualization can be found in the `example_inputs` folder.

.. image:: https://user-images.githubusercontent.com/1488903/61180809-e43cdd80-a61c-11e9-91d7-7d13539c9c16.png
**Clustering of the inverted hexagonal phase with four inner channels connected to a bilayer with a fusion stalk.**

Inside the channels is a fragment of dsDNA. The leaflet clustering was performed using a resolution of 0.5 and hyperesolution turned on. This to allow for the correct clustering of the tight geometry of the channels in coarse grain data (Martini), also force clustering was turned on to have (almost?) every lipid assigned up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61180812-f9b20780-a61c-11e9-838f-f42e54133669.png
**Leaflet clustering of a complex plasmamembrane thether.**

The two leaflets of the plasmamembrane are clearly assigned correctly and depicted as a transparent surface. The cholesterol inside the two leaflets is drawn in VDW spheres and their headgroups have a slightly altering colour. All cholesterol seems to be assigned correctly. Clustering was performed with a 1 nm resolution and forced clustering to assign (all?) the diving cholesterol up to a distance of 2 nm.

.. image:: https://user-images.githubusercontent.com/1488903/61181667-b90cbb00-a629-11e9-9fc0-b2d52e4eaa93.png
**Leaflet clustering of a plasma membrane including multiple proteins.**

The leaflet assignment seemed to have worked correctly, however, we do see some noise in places where we wouldn't expect it at first sight and this behaviour is to be further inspected. For clustering a resolution of 1 nm and forced clustering within 2 nm was used. The protein was used as exclusion to prevent flopping lipids next to proteins to intervene with the leaflet assignment. In total 1.3 millions beads were clustered in less than a minute on a desktop.

Credits
---------

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
