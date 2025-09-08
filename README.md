# NumericalBeamsSU25
Numerical models and analysis for beam vibrations, with the goal of writing a COMAP paper describing the harmonics of a kalimba by adjusting various parameter.

## Mathematica Files
Includes Mathematica files used to verify and test certain analytic solutions of problems such as completeness of FEM. Much of this content was created by my professor, Dr. Allan Struthers.
- **FEMBasisAnalysis.nb** involves analytical analysis of the completeness of different types of FEM bases for both wave and beam problems. Note that much of this file was written by Dr. Struthers.

## Lead In and Basic Boundary Conditions
Includes MATLAB scripts which simulate basic boundary and initial conditions for guitar strings and beams. Solution approaches involve finite difference methods (FDM) and finite element methods (FEM).  
These files include step-by-step descriptions of the solution approach, along with information on specific problems' time-domain solutions and mode shapes over space.

The following files use **FDM** approaches:
- **GuitarStringsFDM.m** contains FDM modeling of a guitar string.
- **BIBIBeamFDM.m** contains FDM modeling of a built in-built in beam.

The following files use **Symbolic FEM** (or analytic FEM) approaches, where FEM techniques are used in conjunction with MATLAB's Symbolic Math Toolbox
- **GuitarStringsSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a guitar string.
- **BIBIBeamSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a built in-built in beam.
- **FFFFBeamSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a free free-free free beam. Its **solution approach is currently incorrect.**
- **BasicKailmbaSymbolicFEM.m** uses symbolic FEM modeling to determine the mode shapes of a basic kalimba.
- **PiecewiseBeamSymbolicFEM.m** uses symbolic FEM modeling to determine if different types of piecewise functions are complete enough to model for beam problems, and to compare to the numeric approach.

The following files use **Numeric FEM** approaches, where FEM is conducted using strictly numerical techniques without leveraging any MATLAB toolboxes:
- **FixedFreeStringFEM.m** uses FEM modeling to determine the mode shapes of a "fixed-free" string. Its **assumptions are incorrect and thus unused.**
- **BIBIBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a built in-built in beam, along with a time-domain solver.
- **BISSBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a built in-simply supported beam, along with a time-domain solver.
- **BIFFBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a built in-free freebeam, along with a time-domain solver.
- **FFFFBeamNumericFEM.m** uses numeric FEM modeling to determine the mode shapes and eigenstructure of a free free-free free beam.
- **PiecewiseBeamNumericFEM.m** uses numeric FEM modeling to determine if different types of piecewise functions are complete enough to model for beam problems.
- **KalimbaNumericFEM.m** uses numeric FEM modeling to determine the mode shapes for a kalimba, along with a time-domain solver.
  
## Functions
Contains MATLAB functions scripts used in other .m files throughout the repository.

## Data Files
Contains prebuilt matrices, colormaps, and other data that can be used in other functions in the repository or to store data that takes a long time to compute. 

## Testbenches
Contains timed tests of different procedures used throughout the repository.
- **FDMTestbench.m** tests different technniques for generating FDM matrices and installing boundary conditions on those matrices. It is **currently outside of the project scope and is unused**.
- **FEMTestbench.m** tests different numerical techniques for generating stiffness matrices for FEM modeling.
- **DITestbench.m** tests different derivation and integration techniques in MATLAB to determine which is most suitable for generating stiffness matrices.
