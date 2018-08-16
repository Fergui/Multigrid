function plot_constr2(X,Y,H,g)
% Call:
% plot_constr2(X,Y,H,g)
%
% Description:
% 2D plot of linear constraints 
%
% Inputs: 
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   H       matrix of linear constraints
%   g       right hand side of the contraints
% Outputs:
%   A 2D plot of the linear constraints
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
% Modification from plot_constr.m by Jan Mandel
%-------------------------------------------------------------------------

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
