# NumericalBeamsSU25
Numerical models and analysis for beam vibrations, with the goal of modeling how to tune beams on a marimba through carving different cross-sectional profiles.

## Lead In and Basic Boundary Conditions
Includes MATLAB scripts which simulate basic boundary and initial conditions for guitar strings and beams. Solution approaches involve finite difference methods (FDM) and finite element methods (FEM).  
These files include step-by-step descriptions of the solution approach, along with information on specific problems' time-domain solutions and mode shapes over space.
- **GuitarStringsFDM.m** contains FDM modeling of a guitar string.
- **GuitarStringFEM.m** contains FEM modeling of a guitar string.
- **BIBIBeamFDM.m** contains FDM modeling of a built in-built in beam.
- **BIBIBeamFEM.m** contains FEM modeling of a built in-built in beam.
- **BasicKailmbaFEM.m** uses FEM modeling to determine the mode shapes of a basic kalimba.
- **FreeFreeBeamFEM.m** uses FEM modeling to determine the mode shapes of a free-free beam.

## Functions
Contains MATLAB functions scripts used in other .m files throughout the repository.

## Data Files
Contains prebuilt matrices, colormaps, and other data that can be used in other functions in the repository or to store data that takes a long time to compute. 

## Testbenches
Contains timed tests of different procedures used throughout the repository.
- **FDMTestbench.m** tests different technniques for generating FDM matrices and installing boundary conditions on those matrices. It is currently outside of the project scope and is unused.
- **IntegrationTestbench.m** tests different numerical integration techniques used for generating stiffness matrices for FEM modeling.
