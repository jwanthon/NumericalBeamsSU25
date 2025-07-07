# NumericalBeamsSU25
Numerical models and analysis for beam vibrations, with the goal of modeling how to tune beams on a marimba through carving different cross-sectional profiles.

## Lead In and Basic Boundary Conditions
Includes MATLAB scripts which simulate basic boundary and initial conditions for guitar strings and beams. Solution approaches involve finite difference methods (FDM) and finite element methods (FEM).  
These files include step-by-step descriptions of the solution approach, along with information on specific problems' time-domain solutions and mode shapes over space.

The following files use **FDM** approaches:
- **GuitarStringsFDM.m** contains FDM modeling of a guitar string.
- **BIBIBeamFDM.m** contains FDM modeling of a built in-built in beam.

The following files use **Symbolic FEM** approaches, where FEM techniques are used in conjunction with MATLAB's Symbolic Math Toolbox
- **GuitarStringsSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a guitar string.
- **BIBIBeamSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a built in-built in beam.
- **FFFFBeamSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a free free-free free beam. Its **solution approach is currently incorrect.**
- **BasicKailmbaSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a basic kalimba.

The following files use **Numeric FEM** approaches, where FEM is conducted using strictly numerical techniques without leveraging any MATLAB toolboxes:
- **FixedFreeStringFEM.m** uses FEM modeling to determine the mode shapes of a "fixed-free" string. Its **assumptions are incorrect and thus unused.**
- **BIBIBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a built in-built in beam.
- **FFFFBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a free free-free free beam.
  
## Functions
Contains MATLAB functions scripts used in other .m files throughout the repository.

## Data Files
Contains prebuilt matrices, colormaps, and other data that can be used in other functions in the repository or to store data that takes a long time to compute. 

## Testbenches
Contains timed tests of different procedures used throughout the repository.
- **FDMTestbench.m** tests different technniques for generating FDM matrices and installing boundary conditions on those matrices. It is **currently outside of the project scope and is unused**.
- **FEMTestbench.m** tests different numerical techniques for generating stiffness matrices for FEM modeling.
- **DITestbench.m** tests different derivation and integration techniques in MATLAB to determine which is most suitable for generating stiffness matrices.
