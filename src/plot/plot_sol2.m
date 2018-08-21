function plot_sol2(X,Y,u,H,g)
% Plot a 3D mesh volume solution u and the constraints Hu=g
%in
%   X,Y   lon and lat matrices
%   u     solution of the fire arrival time
%   H     matrix of interpolation of the perimeters and ignition point
%   g     right hand side of the constraints Hu=g

mesh(X,Y,u), view(2), hold on,
plot_constr_scatter(X,Y,H,g);
end

