# Tests

## Binary search
Test of the binary search routines:
- binsearch_test.m: Test the routine binsearch.m with an easy example.

- binsearch_bounds_test.m: Test the routine binsearch_bounds.m with an easy example.

- binsearch_leq_test.m: Test the routine binsearch_leq.m with an easy example.

## Condensation
Test of the condensation routines:
- condense_test.m: Test the routines addshape.m (to add shape points to the perimeter), interop_bary.m (baryncenters interpolation), and condense.m (condense the matrix H) in a real kml file perimeter using a simulation grid.

## Interpolation
Test of the interpolation routines:
- interop_bary_test.m: Test the routine interop_bary.m with an ideal example.

- interop_bary_real_test.m: Test the routine interop_bary.m with shape files from Las Conchas fire and simulation mesh.

## Initial approximation
Test of the initial approximation routines (Dirichlet problem with FFT sine discretization):
- primalc_test.m: Test the routine primalc.m which solve the saddle point problem with possibly singular symmetric S

			S*u + H'*lambda = f
			H*u             = g

- fdirichlet_fft2_test.m: Test the routine fdirichlet_fft2.m which discretize the minus Laplace operator with dirichlet zero boundary conditions by sine FFT.

- pdirichlet_constr_test.m: Test the routines interop_bary.m and apdirichlet_constr.m in an ideal case of two concentric perimeters. The pdirichlet_constr.m routine solves the initial approximation u formulated where L = d^2/dx^2 + d^2/dy^2 is the Laplace operator as

			(-L)^a * u = 0 on rectangle 
        	u          = bv on the boundary, subject to   
			H*u        = g

## Setup
Test the setup routines:
- setup_real_test.m: Test the routine setup_real.m. An input file in data/in.mat where the preprocessing part is done using preproc.m routine is required. It saves a Matlab file called setup.mat with the setup of the case for the optimization.

## Multigrid
Test the Multigrid routines:
- multigrid_ideal_test.m: Test the routines setup_ideal.m and multigrid_process.m. The routine setup_ideal.m setups the ideal case from mesh dimensions and spacing. The routine multigrid_process.m runs the projected Multigrid descent method.

- multigrid_ideal_test_q.m: Test the same routines than the previous multigrid_ideal_test.m, but using different values for p-norm in the objective function.

- multigrid_file_test.m: Test the routines setup_file.m and multigrid_process.m. The routine setup_file.m setups the ideal case from WRF-SFIRE simulation from a wrfout file in data/wrfout and the simulated frame of the perimeter. The routine multigrid_process.m runs the projected Multigrid descent method.

- multigrid_real_test.m: Test the routines preproc.m, setup_real.m and multigrid_process.m. The routine preproc.m preprocess the data from the wrfout simulation files in the path sim_path, the KML files of the perimeters in the path per_path, the ignition time in itime, the ignition point coordinates in ilon and ilat, and using or not using dynamic rate of spread in dyn. This preproc.m generates an in.mat file which is used for the routine setup_real.m for setup the real case. This setup_real.m generates a file called setup.mat. Finally, the routine multigrid_process.m runs the projected Multigrid descent method. 

To run a real example, the next steps are required changing the script multigrid_real_test.m:
1) Remove all the files inside data (in.mat and setup.mat) from previous executions (because if it exists data inside it is going to load this data directly and not computes again the same thing).
2) Change the paths to the simulation files and KML files.
3) Change the ignition time (using the same format).
4) Change the latitude and longitude coordinates of the ignition point.
5) Decide if to use or not a dynamic rate of spread. Warning: using dynamic rate of spread, a huge memory is required.
6) Run multigrid_real_test.m in Matlab.
7) Wait for the results in a Matlab file called out.mat!
