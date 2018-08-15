function plot_constr(X,Y,H,g,m,bounds)
%plot_constr(X,Y,H,g,m,bounds)
% 3D plot of linear constraints
%in 
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   H       matrix of linear constraints
%   g       right hand side of the contraints
%   m       type of plot (example: 'r.' for red dots)
%   bounds  bound where plotting the data
%out
%   A 3D plot of the linear constraints

[nx,ny]=size(X);
if nx==1,
    nx=X;
    ny=Y;
    [X,Y]=ndgrid(1:nx,1:ny);
end
% constraint nodes
xx=H*X(:); 
yy=H*Y(:);
plot3(xx,yy,g,m);
ix=find(xx >= bounds(1) & xx <= bounds(2) & yy >= bounds(3) & yy <= bounds(4));
plot3(xx(ix),yy(ix),g(ix),m);
