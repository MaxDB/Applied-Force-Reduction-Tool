# Applied Force Reduction tool
- AFR tool is designed to produce reduced-order models of geometrically nonlinear systems.
- Whilst primarily aimed at finite element models, analytically defined models can also be analysed. 
- The method is described in detail in '[Evaluating applied force reduced-order models: verification](https://doi.org/10.1007/s11071-026-12899-6)' and '[Evaluating applied force reduced-order models: validation](https://doi.org/10.1007/s11071-026-12844-7)'.

## Installation and use
- Download and extract or clone the repository.
- The directory containing the repository can be renamed and placed anywhere.
- Open the project file: `AFR_tool.prj` in MATLAB. This will temporarily add all the required directories to the path.
- The project file needs to be open for the tool to run, but once open, it can be used from any directory.

## Getting started
- Several 'Example.m' live scripts provide walkthroughs for reproducing published results.

## Required and recommended dependencies
### MATLAB toolboxes
- Optimization Toolbox: required.
- Symbolic Math Toolbox: required for analytically defined models.
- Parallel Computing Toolbox: recommended for increased performance.

### External software
- Continuation Core (COCO): https://sourceforge.net/projects/cocotools/.
- Abaqus (the free learning edition is also supported https://www.3ds.com/edu/education/students/solutions/abaqus-le).  

## Acknowledgements
- Colour palettes for plots are taken from colorbrewer2: https://colorbrewer2.org/.
