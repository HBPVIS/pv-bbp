-------------------------------------------------------------------------
WHAT IS IT ?
-------------------------------------------------------------------------
Paraview plugin to read the BBP Meshes using the BBP-SDK Mesh Parser.

-------------------------------------------------------------------------
HOW TO BUILD IT ?
-------------------------------------------------------------------------
- edit the proj/cmake/configure to set the path to the Paraview BUILD directory
- compile and link the plugin as a shared library 

cd proj/cmake
./configure
make

-------------------------------------------------------------------------
HOW TO USE IT ?
-------------------------------------------------------------------------
- start paraview
- go to Tools > Plugin Manager
- press "Load New" and load the plugin produced before (.so file in lib directory) 
- open a mesh file (exemple in data directory)
- display the "Object Inspector" and click on "Apply" to see the neuron.

-------------------------------------------------------------------------
EXTRA DOCUMENTATION
-------------------------------------------------------------------------
https://bbpteam.epfl.ch/wiki/index.php/BBP_Paraview

