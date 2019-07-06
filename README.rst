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

This software has been developed to allow for a higher selection syntax than atom and or residue index. MDVoxelClustering combines nieghbour clustering with a voxel mask of set resolution (default 1 nm). For now the software has mainly been tested on CG Martini lipid assemblies with or witout embedded proteins, though it probably works even better on atomistic simulations due to the higher amount of atoms that can be binned. 

* Free software: ISC license
* Documentation: https://mdvoxelclustering.readthedocs.org.

Features
--------
* This is an alpha build and should be used with extreme care!
* Voxel based neighbour clustering
* Fast contour clustering
* Compatible with most MD file formats due to its tight link to MDAnalysis
* Consistent clustering over time
* Membrane leaflet recognition of lipids of most topologies
  - Bilayers
  - Vesicles
  - Inverted hexagonal phase
  - Membrane thethers
  - Complex lipids formulations including cholesterol


Credits
---------

Tools used in rendering this package:

*  Cookiecutter_
*  `cookiecutter-pypackage`_

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
