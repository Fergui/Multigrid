function plot_sol2(X,Y,u,H,g)
% Call:
% plot_sol2(X,Y,u,H,g)
%
% Description:
% Plot a 2D mesh volume solution u and the constraints Hu=g
%
% Inputs:
%   X,Y   lon and lat matrices
%   u     solution of the fire arrival time
%   H     matrix of interpolation of the perimeters and ignition point
%   g     right hand side of the constraints Hu=g
% Outputs:
%   A 2D plot of the solution u with the constraints Hu=g
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

mesh(X,Y,u), view(2), hold on,
plot_constr_scatter(X,Y,H,g);
end

