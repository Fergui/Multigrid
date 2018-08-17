function [un,Jop]=locallinesearch(u,node,phi,Jop,p)
% Call:
% um=locallinesearch(node,phi,p)
%
% Description:
% Computes the exact line search in the patch 
%     j+s
%i-s       i+s
%     j-s
% using a tend base function for the direction
%
% Inputs:
%   u      fire arrival time
%   node   array containing
%       i    x node index of the center of the patch
%       j    y node index of the center of the patch
%       s    spacing of the patch
%   phi    coarse base function
%   Jop    array of objective function values
%   p      structure
% Outputs:
%   un     final result of a local linesearch
%   Jop    adding the new objective function value
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

[m,n]=size(u);
i=node(1);
j=node(2);
s=node(3);
% Defining the computing bounds: [li,ui] and
% [lj,uj] and the result bounds: [lir,uir] and
% [ljr,ujr]. For computing the solution in the
% boundaries of the patch.
li=max(i-s-1,1);
lj=max(j-s-1,1);
ui=min(i+s+1,m);
uj=min(j+s+1,n);
mms=size(li:ui,2);
nns=size(lj:uj,2);
lir=i-s;
ljr=j-s;
uir=i+s;
ujr=j+s;
si=1+lir-li:mms+uir-ui;
sj=1+ljr-lj:nns+ujr-uj;
p.bi=si;
p.bj=sj;
% Columns of H in the patch
cols=sum(combvec((li:ui),((lj:uj)-1)*m),1);
% H matrix in the specific patch
H=p.H;
H=H(:,cols);
ind=(sum(H,2)<eps);
H(ind,:)=[];
g=p.g;
g(ind)=[];
% Compute the part of u where we are going to compute the optimization
um=u(li:ui,lj:uj);
% Compute the part of R where we are going to compute the optimization
R=p.R(li:ui,lj:uj);
% Declaring the direction
dir=zeros(mms,nns);
% Defining the direction from phi
dir(1-li+lir:end-ui+uir,1-lj+ljr:end-uj+ujr)=phi;
% Project the direction into Hd=0
d=projectdir(dir,H);
% Looking to the negative and positive direction
for sig=-1:2:1
    dir=sig*max(0,d);
    % Line search in the patch
    [Jmin,~,um,~]=linesearch(um,@cJ,R,dir,p);
    % Saving the results of the line search
    if Jmin<Jop(end)
        Jop=[Jop;Jmin];
    else
        Jop=[Jop;Jop(end)];
    end
end
un=u;
un(i-s:i+s,j-s:j+s)=um(si,sj);
end

