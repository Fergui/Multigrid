function plot_constr_scatter(X,Y,H,g)
% Call:
% plot_constr(X,Y,H,g)
%
% Description:
% 3D plot of linear constraints
%
% Inputs: 
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   H       matrix of linear constraints
%   g       right hand side of the contraints
% Outputs:
%   A 3D plot of the linear constraints
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
% Modified from plot_constr.m by Jan Mandel
%-------------------------------------------------------------------------

uu=unique(g);
[nx,ny]=size(X);
if nx==1,
    nx=X;
    ny=Y;
    [X,Y]=ndgrid(1:nx,1:ny);
end
% constraint nodes
xx=H*X(:); 
yy=H*Y(:);
zz=g;
col=[255,0,0;255,165,0;0,100,255]/255;
c=ones(length(zz),3);
for k=1:length(uu)
    c(zz==uu(k),:)=repmat(col(k,:),length(find(zz==uu(k))),1);
end
scatter3(xx,yy,zz,5,c,'*');
