Requirements:

C++ compiler supporting OpenMP (works with g++)
C++ version > 17, compile with -std==c++17
if C++ version is not > 17, possibly would cause Eigen Alignment Issue, the code has been modified to consider it
but it is recommended to use c++17 or above to circumvent it see:https://eigen.tuxfamily.org/dox/group__DenseMatrixManipulation__Alignement.html
MPI libraries (works with OpenMPI and MPICH)
Eigen library (http://eigen.tuxfamily.org/index.php?title=Main_Page)

Header files:

definitions.h - includes, etc
Grain3d.h - the Grain3d class
Levelset3d.h - Levelset3d class, each Grain3d has a Levelset3d as a member.  All level sets are in 1-dimensional vectors that go along x first, then y, then z so the access is lset(x,y,z) = lset(z*ydim*xdim + y*xdim + x).  This information is only really important in making sure the input level set (from an input file) is in the format that the Levelset3d requires.
readParams.h - reads input files
Wall3d.h - Wall3d class, infinite plane walls

Currently, the MainTest.cpp can be compiled with:
mpiCC -std=c++17 -I /path/to/eigen3 MainTest.cpp -O3 -o maintest
