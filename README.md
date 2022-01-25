# MDVoxelSegmentation
Using neighbor segmentation in voxelspace for fast and consistent spatial and temporal segmentation.

This software has been developed to allow for a higher selection syntax than atom and or residue index, such as abstract complex particles (e.g. lipid monolayers). MDVoxelSegmentation combines multiple layers of neighbor segmentation with a voxel mask, and makes tracking segments over time possible at high quality.

This code was published in this open access [https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c00446](article).

* Open software: Apache 2 license

## Features
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

