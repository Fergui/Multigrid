function plot_constr2(X,Y,H,g)
%plot_constr2(X,Y,H,g)
% 2D plot of linear constraints 
%in 
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   H       matrix of linear constraints
%   g       right hand side of the contraints
%out
%   A 2D plot of the linear constraints

uu=unique(g);
[nx,ny]=size(X);
if nx==1,
    nx=X;
    ny=Y;
    [X,Y]=ndgrid(1:nx,1:ny);
end
% constraint nodes
col=[255,0,0;255,165,0;0,100,255]/255;
xx=H*X(:); 
yy=H*Y(:);
for l=1:length(uu)
    k=find(g==uu(l));
    plot(xx(k),yy(k),'.','Color',col(l,:)), hold on
end