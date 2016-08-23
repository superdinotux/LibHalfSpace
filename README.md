# LibHalfSpace
code for the paper "LibHalfSpace: a C++ object-oriented library to study deformation and stress in elastic half-spaces"


LibHalfSpace is a coherent programming interface containing a set of tools designed to make easier and faster the study of processes in an elastic half-space. 
LibHalfSpace is presented in the form of an object-oriented library.
A set of well known and frequently used source models (Mogi source, penny shaped horizontal crack, inflating spheroid, Okada rectangular dislocation, etc.) are implemented to describe the potential usage and the versatility of the library.
The common interface given to library tools enables us to switch easily among the effects produced by different deformation sources that can be monitored at the free surface.
Furthermore, the library also offers an interface which simplifies the creation of new source models exploiting the features of object-oriented programming (OOP).
These source models can be built as distributions of rectangular boundary elements.
In order to better explain how new models can be deployed some examples are included in the library.

The project presents the following directories:
- Documentation, the directory which contains the library documentation
- Examples, a directory which contains the sources of the examples described in the paper
- SRC, the directory with the C++ source code
- SRC_fortran, this directory contains the original dc3d.f Okada code
