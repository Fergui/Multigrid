function plot_sol_mesh(X,Y,u,H,g)
% Call:
% plot_sol_mesh(u,H,g)
%
% Description:
% Plot the fire arrival time and the constraints in a 3D contour and mesh
%
% Inputs:
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   u   fire arrival time over the domain
%   H   matrix of the constraints Hu=g
%   g   right hand side of the constraints Hu=g
% Outputs:
%   Mesh 3D plot of the solution u and the constraints Hu=g
% 
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

[nx,ny]=size(X);
if nx==1,
    nx=X;
    ny=Y;
    [X,Y]=ndgrid(1:nx,1:ny);
end

mesh(X,Y,u), hold on
plot_constr_scatter(X,Y,H,g);
end
