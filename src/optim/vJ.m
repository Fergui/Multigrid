function [C,r] = vJ(u,R,p)
% Call:
% C = vJ(u,R,p) 
%
% Description:
% Evaluates J(u,R,p) and returns the contribution matrix
%
% Inputs:
%   u         level set function
%   R         rate of spread, on same nodes as u
%   p         parameter structure with fields:  
%                dx,dy   mesh spacing
%                select  handle to upwinding function
%                ofunc   objective function f(x,y) comparing 
%                        x=||grad u||^2 and y=R^2 where xy=1 
% Outputs:
%   C         contribution matrix
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

[m,n]=size(u);
diffLx=zeros(m,n);
diffRx=zeros(m,n);
diffLy=zeros(m,n);
diffRy=zeros(m,n);
diffLx(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(1:end-2,2:end-1))/p.dx;
diffRx(2:end-1,2:end-1)=(u(3:end,2:end-1)-u(2:end-1,2:end-1))/p.dx;
diffLy(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(2:end-1,1:end-2))/p.dy;
diffRy(2:end-1,2:end-1)=(u(2:end-1,3:end)-u(2:end-1,2:end-1))/p.dy;
diffx=p.select(diffLx,diffRx,p.dx);
diffy=p.select(diffLy,diffRy,p.dy);
sngrad=diffx.^2+diffy.^2;
ros=R.^2;
C=p.ofunc(sngrad,ros);
C(~p.vmask)=0;
C=C(2:end-1,2:end-1);
r=(sum(C(:).^p.q)*p.dx*p.dy)^(1/p.q);

end
