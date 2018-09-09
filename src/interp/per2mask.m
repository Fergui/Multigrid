function [M,u]=per2mask(H,g,m,n,dx,dy)
%M=per2mask(H)
% Define a mask from perimeter defined as constraints Hu=g
%in 
%   H      matrix of constraints
%   g      right hand side of the constraints
%   m,n    x and y dimensions
%   dx,dy  x and y resolution
%out
%   M      mask matrix (1's inside and 0's outside)
%   u      solution of Dirichlet boundary conditions using constraints

uu=unique(g);
gs=[m,n];
h=[dx,dy];
a=1.4;
bv=uu(end)*10;
u=pdirichlet_constr(gs,h,bv,H,g,a);
M=(u<=uu(end));

end