function [M,u]=per2mask(H,g,m,n,dx,dy)
% Call:
% M=per2mask(H)
%
% Description:
% Define a mask from perimeter defined as constraints Hu=g
%
% Inputs: 
%   H      matrix of constraints
%   g      right hand side of the constraints
%   m,n    x and y dimensions
%   dx,dy  x and y resolution
% Outputs:
%   M      mask matrix (1's inside and 0's outside)
%   u      solution of Dirichlet boundary conditions using constraints
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

uu=unique(g);
gs=[m,n];
h=[dx,dy];
a=1.4;
bv=uu(end)*10;
u=pdirichlet_constr(gs,h,bv,H,g,a);
M=(u<=uu(end));

end
