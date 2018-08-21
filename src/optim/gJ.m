function [varargout] = gJ(u,p)
% r = gradJ(u,p) 
%           Evaluates J(u,p)
% [gJp,gJn] = gradJ(u,p)
%           Evaluates the gradient of J from the positive and negative
%           side
% [gJp,gJn,sngrad] = gradJ(u,p)
%           Evaluates the gradient of J from the positive and negative side
%           and ||grad u||^2
%inputs:
% u         level set function
% p         parameter structure with fields:
%                vmask   matrix, true where the values of level set 
%                        function can change  
%                R       rate of spread, on same nodes as u
%                dx,dy   mesh spacing
%                select  handle to upwinding function
%                ofunc   objective function f(x,y) comparing 
%                        x=||grad u||^2 and y=R^2 where xy=1 
%                dfdG    partial derivative of objective function respect 
%                        to ||grad u||^2
%                dfdR    partial derivative of objective function respect 
%                        to R(u)^2
%outputs:
% r         Value of J using u and p.
% gJn       4-D matrix, d(J(u))/du in positive direction
% gJp       4-D matrix, d(J(u))/du in negative direction
% sngrad    2-D matrix, ||grad u||^2

[m,n]=size(u);
diffLx=zeros(m,n);
diffRx=zeros(m,n);
diffLy=zeros(m,n);
diffRy=zeros(m,n);
diffLx(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(1:end-2,2:end-1))/p.dx;
diffRx(2:end-1,2:end-1)=(u(3:end,2:end-1)-u(2:end-1,2:end-1))/p.dx;
diffLy(2:end-1,2:end-1)=(u(2:end-1,2:end-1)-u(2:end-1,1:end-2))/p.dy;
diffRy(2:end-1,2:end-1)=(u(2:end-1,3:end)-u(2:end-1,2:end-1))/p.dy;
[diffx,diffGxdu]=p.select(diffLx,diffRx,p.dx);
[diffy,diffGydu]=p.select(diffLy,diffRy,p.dy);
sngrad=diffx.^2+diffy.^2;
ros=p.R.^2;
c=p.ofunc(sngrad,ros);
c(~p.vmask)=0;
%figure(10), mesh(c');
r=(sum(c(:).^p.q)*p.dx*p.dy)^(1/p.q);
varargout{1}=r;
if nargout >= 2
    dfdG=zeros(m,n);
    dfdR=zeros(m,n);
    dfdG(2:end-1,2:end-1)=p.dfdG(sngrad(2:end-1,2:end-1),ros(2:end-1,2:end-1));
    dfdR(2:end-1,2:end-1)=p.dfdR(sngrad(2:end-1,2:end-1),ros(2:end-1,2:end-1));
    dfdG(~p.vmask)=0;
    dfdR(~p.vmask)=0;
    sn=zeros(m,n);
    sp=zeros(m,n);
    Kx=abs(c(2:end-1,2:end-1)).^(p.q-1).*sign(c(2:end-1,2:end-1)).*dfdG(2:end-1,2:end-1).*(2*diffx(2:end-1,2:end-1));
    Ky=abs(c(2:end-1,2:end-1)).^(p.q-1).*sign(c(2:end-1,2:end-1)).*dfdG(2:end-1,2:end-1).*(2*diffy(2:end-1,2:end-1));
    sn(2:end-1,2:end-1)=sn(2:end-1,2:end-1)+Kx.*diffGxdu.cn(2:end-1,2:end-1)+Ky.*diffGydu.cn(2:end-1,2:end-1);
    sn(1:end-2,2:end-1)=sn(1:end-2,2:end-1)+Kx.*diffGxdu.ln(2:end-1,2:end-1);
    sn(3:end,2:end-1)=sn(3:end,2:end-1)+Kx.*diffGxdu.rn(2:end-1,2:end-1);
    sn(2:end-1,1:end-2)=sn(2:end-1,1:end-2)+Ky.*diffGydu.ln(2:end-1,2:end-1);
    sn(2:end-1,3:end)=sn(2:end-1,3:end)+Ky.*diffGydu.rn(2:end-1,2:end-1);
    sp(2:end-1,2:end-1)=sp(2:end-1,2:end-1)+Kx.*diffGxdu.cp(2:end-1,2:end-1)+Ky.*diffGydu.cp(2:end-1,2:end-1);
    sp(1:end-2,2:end-1)=sp(1:end-2,2:end-1)+Kx.*diffGxdu.lp(2:end-1,2:end-1);
    sp(3:end,2:end-1)=sp(3:end,2:end-1)+Kx.*diffGxdu.rp(2:end-1,2:end-1);
    sp(2:end-1,1:end-2)=sp(2:end-1,1:end-2)+Ky.*diffGydu.lp(2:end-1,2:end-1);
    sp(2:end-1,3:end)=sp(2:end-1,3:end)+Ky.*diffGydu.rp(2:end-1,2:end-1);
    gJn=r^(1-p.q)*sn;
    gJp=r^(1-p.q)*sp;
    varargout{1}=gJn;
    varargout{2}=gJp;
    if nargout>=3
        varargout{3}=sngrad;
    end
end
end
