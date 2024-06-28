# Applied Force Reduction Tool
- AFR tool is designed to produce reduced-order models of geometrically nonlinear systems.
- Whilst primarily aimed at finite element models, analytically defined models can also be analysed. 
- The method is described in detail in [to be published].

## Installation and Use
- Download and extract or clone the repository.
- The directory containing the repository can be renamed and placed anywhere.
- Open the project file: `AFR_tool.prj` in MATLAB. This will temporarily add all the required directories to the path.
- The project file needs to be open for the tool to run, but once open, can be used from any directory.

## Required and Recommended Dependencies
### MATLAB Toolboxes
- Optimization Toolbox: required
- Symbolic Toolbox: required for analytically defined models
- Parallel Toolbox: recommended for increased performance

One of each of the following dynamic and FE solvers is also required
### Supported Dynamic Solvers
- Continuation Core: 
  - Free and open source
  - https://sourceforge.net/projects/cocotools/

### Supported FE Solvers
- Abaqus: 
  - Not free
  - https://www.3ds.com/products/simulia
- Abaqus Learner Edition:
  - Free
  - Limited to 1,000 nodes
  - https://www.3ds.com/edu/education/students/solutions/abaqus-le

## Acknowledgements
- Colour palettes for plots are taken from colorbrewer2: https://colorbrewer2.org/
