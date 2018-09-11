# Multigrid

Matlab framework to compute projected preconditioned gradient descent method and projected multigrid descent method. Experiments optimizing the residues of the Eikonal equation for fire arrival time in a coupled atmosphere-fire data (from WRF-SFIRE simulations). File startup.m includes into the path all the source routines. 

## Source

### fft
Routines to compute sine FFT for the minus Laplacian operator with Dirichlet boundary conditions.
### fire
Routines to compute fire dynamics as rate of spread, fuels, and interpolation of wind components and fuel moisture.
### interp
Routines to interpolate and condense operator H in the linear system of constraints Hu=g.
### netcdf
Routines to read data from Netcdf files (WRF-SFIRE simulation output files).
### optim
Routines in the optimization process: computation of the objective function, binary search, linesearch, projected preconditioned gradient descent method and projected multigrid descent method.
### plot
Routines to plot results.
### setup
Routines to setup cases: ideal, ideal from WRF-SFIRE simulation and real cases.
### utils
Other routines as: abstract KML files information to matlab, record the graphics in a gif file, etc.

## Tests

To test the functions:
1) Start Matlab inside Multigrid folder, so inside Multigrid folder type in terminal

      \$ matlab -nodesktop -nodisplay

2) If you started Matlab in another folder, run inside Matlab 

      \>\> startup

3) Run any of the tests in the tests folder. For example:

      \>\> cd tests

      \>\> pdirichlet_constr_test
