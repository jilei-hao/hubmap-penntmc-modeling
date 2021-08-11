# hubmap-penntmc-modeling

This repository includes 3D modeling code developed for the Penn Center for Multi-scale Molecular Mapping of the Female Reproductive System in the HuBMAP Consortium.

Tools:
OvaryModelGenerator is source code for creating a multi-label 3D multi-label image of an ovary, subdivided into long-axis and rotational subregions. 
To build the executable from source code, you need:
- CMake, the multi-platform build environment from http://cmake.org
- VTK, the Visualization Toolkit from http://vtk.org

The building procedure is:
- Install cmake
- Build VTK using cmake
- Build OvaryModelGenerator using cmake
