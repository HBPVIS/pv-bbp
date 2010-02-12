-------------------------------------------------------------------------
WHAT IS IT ?
-------------------------------------------------------------------------
Paraview plugin to read the BBP HDF5 Morphologies using the BBP-SDK Morphology Parser.

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
- open a morphology file (exemple in data directory)
- display the "Object Inspector" and click on "Apply" to see the neuron.

Notes : You may have some troubles with the HDF5 library while loading the file. If this 
is the case, you should recompile Paraview using the hdf5 library from the system and not 
the one from VTK. To do that, change the option PARAVIEW_USE_SYSTEM_HDF5 to ON in the
configuration of the Cmake before the compilation of Paraview.

-------------------------------------------------------------------------
EXTRA DOCUMENTATION
-------------------------------------------------------------------------
https://bbpteam.epfl.ch/wiki/index.php/BBP_Paraview

