function plot_sol(X,Y,u,H,g)
% Call:
% plot_sol(X,Y,u,H,g)
%
% Description:
% Plot the fire arrival time and the constraints in a 3D contour
%
% Inputs:
%   X,Y   lon and lat matrices
%   u   fire arrival time over the domain
%   H   matrix of the constraints Hu=g
%   g   right hand side of the constraints Hu=g
% Outputs:
%   3D contour plot of the solution u with the constraints Hu=g.
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

uu=unique(g);
contour3(X,Y,u,[uu(1):(uu(end)-uu(1))/10:uu(end)],'k'), hold on
plot_constr_scatter(X,Y,H,g);
end

