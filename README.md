# MDVoxelClustering 
#### An attempt of making data analysis easier and quicker for many MD people

## MAIN GOAL
To do analysis of a MD trajectory on a more abstract level. Instead of having to talk about specific atoms and their corresponding position and residue members. We often would like to talk about properties on a much more abstract level, such as the diffusion of a certain aggregate of particles, or talk about the leaflets of lipid assemblies. Current methods for leaflet identification often imply geometrical arguments  on every residue, which can quickly become very expensive. Here we present a mixture of voxel and graph based selection methods to obtain information about dynamic clusters at a (hopefully) realtime speed.

### Clustering
#### CONTOURS
Contours are meant for analysis on the behaviour of clusters in a cheaper dimensionality. They could be enough to uniquely identify a cluster e.g. the number of particles present and the voxel space contour size (not exactly the real size, but often close enough).

#### VOLUMES
Volumes give you a more robust manner of selection which will for sure include all the particles you need for high resolution analysis. Combining layers of volumes and contours often can dramatically reduce the degrees of freedom in your data set, without loosing any particles on the way. Hopefully we will demonstrate this for the automatic leaflet detection in a crowded membrane space (all forms of lipid cluster states, such as: adhesed, semi-fused, fused and seperated). Including densely populated membranes with proteins. 

