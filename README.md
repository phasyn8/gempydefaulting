"# gempydefaulting" 

IMPLICIT Geological modeling framework GemPy allows geologist and modelers to make and evaluate structural geologic models for use
in exploration, research, and subsurface investigations. Experienced users will know that the software is powerful in its
inplementation and ability to create useful models. A limitaiton for this software is the way that it handles the interface between 
faults and structural discontinuities. As the 3D visualizations of the models are distillations of a scalar field they are subject 
to artifaces that arise from the structures grids that make up the modeling volume. These lead to stair-stepping artifacts that
occur partiucularly often when formations are displaced by faults.

This module is an attempt to patch this behaviour, only for visualization, by eliminating these artifacts and patching the gaps
so as to firstly clean up the visualization as well as conform to a more geologically realistic subsurface model. The goal is to
provide a simple function that will take a computed GemPy model and pass through this filter to create a 3D visualization. 
In addition to Gempy, this package also utilizes PyVista for visualization and GemGIS for surface generation.


Please direct comments and features requests this git repo.

phasyn8 - 24/04/23
