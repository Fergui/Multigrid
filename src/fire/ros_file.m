function ros=ros_file(u,ign,per,p)
% Call:
% ros=ros_file(u,ign,per,p)
%
% Example: 
% ros=ros_file(u,ign,per,p)
%
% Description:
% Computes rate of spread from WRF-SFIRE output file in structures
%
% Inputs:
%   u      Fire arrival time
%   ign    Structure from the ignition time, with:
%               dzdxf   x component of the slope
%               dzdyf   y component of the slope
%   per    Structure from the first perimeter, with:
%               uf      x component of the wind
%               vf      y component of the wind
%               fmc_g   fuel moisture in the first perimeter
%   p      Static structure, with:
%               dx, dy       fire mesh spacing
%               nfuelcat     matrix of vegetation types 
% Outputs:
%   ros    Rate of spread computed from WRF-SFIRE output file
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Take input data in s structure
dx=p.dx;
dy=p.dy;
%% First declarations
[m,n]=size(u);
diffRx=zeros(m,n);
diffLx=zeros(m,n);
diffRy=zeros(m,n);
diffLy=zeros(m,n);
%% Normal direction nvx and nvy 
diffRx(2:end-1,2:end-1)=(u(3:end,2:end-1)-u(2:end-1,2:end-1))/dx;
diffLx(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(1:end-2,2:end-1))/dx;
diffRy(2:end-1,2:end-1)=(u(2:end-1,3:end)-u(2:end-1,2:end-1))/dy;
diffLy(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(2:end-1,1:end-2))/dy;
diffCx=(diffRx+diffLx)/2;
diffCy=(diffRy+diffLy)/2;
scale=sqrt(diffCx.*diffCx+diffCy.*diffCy+eps);
nvx=diffCx./scale;
nvy=diffCy./scale;
%% Wind speed in spread direction
speed=per.uf.*nvx+per.vf.*nvy;
%% Slope in spread direction
tanphi=ign.dzdxf.*nvx+ign.dzdyf.*nvy;
%% Rate of Spread
ros=fire_ros(p.nfuelcat,speed,tanphi,per.fmc_g);
end