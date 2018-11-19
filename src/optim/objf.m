function [C,r]=objf(u,R,p)
% Call:
% [C,r] = objf(u,R,p) 
%
% Description:
% Evaluates
%   integral (sqrt(|u_xx|^2 +  |u_yy|^2 + |u_xy|^2 + |u_yx|^2)),
% and returns the contribution matrix C and the objective function value r
%
% Inputs:
%   u         level set function
%   p         parameter structure with fields:  
%                dx,dy   mesh spacing
%                select  handle to upwinding function
%
% Outputs:
%   C         contribution matrix
%   r         objective function value
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-11-15
%-------------------------------------------------------------------------

[m,n]=size(u);
dx=p.dx;
dy=p.dy;
u_xx=zeros(m,n);
u_yy=zeros(m,n);
u_xy=zeros(m,n);
u_yx=zeros(m,n);
u_xx(2:end-1,2:end-1)=( (u(3:end,2:end-1)-u(2:end-1,2:end-1))/dx - (u(2:end-1,2:end-1)-u(1:end-2,2:end-1))/dx )/dx;
u_yy(2:end-1,2:end-1)=( (u(2:end-1,3:end)-u(2:end-1,2:end-1))/dy - (u(2:end-1,2:end-1)-u(2:end-1,1:end-2))/dy )/dy;
u_xy(2:end-1,2:end-1)=( (u(3:end,3:end)-u(1:end-2,3:end))/(2*dx) - (u(3:end,1:end-2)-u(1:end-2,1:end-2))/(2*dx) )/(2*dy);
u_yx(2:end-1,2:end-1)=( (u(3:end,3:end)-u(3:end,1:end-2))/(2*dy) - (u(1:end-2,3:end)-u(1:end-2,1:end-2))/(2*dy) )/(2*dx);
C=sqrt(u_xx.*u_xx+u_yy.*u_yy+u_xy.*u_xy+u_yx.*u_yx);
C=C(2:end-1,2:end-1);
r=sum(C(:));

end