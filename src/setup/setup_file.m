function [u,p] = setup_file(file,perF)
%[u,p] = setup_file(file,perT)
%example: [u,p] = setup_file('wrfout',50);
% Set up u and p for the line search in the case of ideal simulation of
% WRF-SFIRE
%input 
%   file     NetCDF WRF-SFIRE output file
%   perF     wrfout second perimeter frame
%output
%   u        level set initialized using an approximation
%   p        p structure with:
%               perF        wrfout perimeter frame
%               lastF       wrfout last frame
%               ignF        wrfout ignition frame
%               dx, dy      fire mesh spacing
%               per1_mask   matrix, mask of the first perimeter
%               per12_mask  matrix, mask of the second perimeter line
%               per2_mask   matrix, mask of the second perimeter
%               per1_time   fire arrival time of the first perimeter
%               per2_time   fire arrival time of the second perimeter
%               mask        matrix, true where the values of level set function have to change
%               vmask       matrix, true where the values of level set function can change
%               bc          boundary conditions = fixed values of the level set function where mask is false
%               nfuelcat    structure with all the fuel data necessary to
%                           compute ROS
%               R           matrix, rate of spread, on same nodes as u
%               ofunc       matlab function, objective function comparing 
%                           x=||grad u||^2 and y=R^2 such that xy=1
%               dfdG        matlab function, partial derivative of ofunc to
%                           respect to x=||grad u||^2
%               dfdG        matlab function, partial derivative of ofunc to
%                           respect to y=R^2
%               q           q norm of the computation of J
%               h           the stepsize to compute the gradient
%               stepsize    step size for minimization
%               numsteps    number of steps to try
%               max_iter    maximum number of iterations
%               select      handle to upwinding function
%               max_step    max size of the steps to search
%               nmesh       number of mesh points in each search
%               max_depth   max number of searchs
%               min_depth   min number of searchs
%               umax        array, maximal value of u
%               umin        array, minimal value of u

%% Take data from files
s=nc2struct(file,{'Times','FIRE_AREA'},{});
times=char(s.times)';
stri=times(1,:);
ns=size(times,1);
simt=zeros(ns,1);
for k=1:ns
    strf=times(k,:);
    simt(k)=str2times(stri,strf)*60;
end
p.perF=perF;
p.lastF=ns;
for t=1:p.lastF
    if max(max(s.fire_area(:,:,t))) > 0
        p.ignF=t;
        break;
    end
end
ignS=nc2struct(file,{'FXLONG','FXLAT','XLONG','XLAT','FIRE_AREA','NFUEL_CAT','DZDXF','DZDYF','UF','VF','FMC_G'},{'DX','DY'},p.ignF);
perS=nc2struct(file,{'FIRE_AREA','TIGN_G','ROS'},{},p.perF);
perlS=nc2struct(file,{'FIRE_AREA','TIGN_G','ROS'},{},p.perF+20);
dynS=nc2struct(file,{'UF','VF','FMC_G'},{});
dynS.simt=simt;
ts=simt(p.ignF);
te=simt(p.perF);
tel=simt(p.perF+20);
%% Number of points and spacing
[m,n]=size(ignS.fire_area);
sr_x=round(size(ignS.fxlong,1)/size(ignS.xlong,1),1);
sr_y=round(size(ignS.fxlat,1)/size(ignS.xlat,1),1);
p.dx=ignS.dx/sr_x;
p.dy=ignS.dy/sr_y;
%% Perimeters
per1=logical(((ignS.fire_area>0).*(ignS.fire_area<1)));
[xx,yy,~]=find(per1);
xx=floor(mean(xx));
yy=floor(mean(yy));
per1=zeros(m,n);
per1(xx,yy)=1;
per12=logical(((perS.fire_area>0).*(perS.fire_area<1)));
per2=logical(perS.fire_area<1);
per3=logical(perlS.fire_area<1);
%% Masks
p.per1_mask=logical(per1);
p.per12_mask=per12;
p.per2_mask=per2;
p.per3_mask=per3;
p.per1_time=ts;
p.per2_time=te;
p.per3_time=tel;
p=setup_masks(p);
%% Shape points
X=ignS.fxlong;
Y=ignS.fxlat;
xi=xx*p.dx;
yi=yy*p.dy;
xq=X(per12);
yq=Y(per12);
%% Matrix H and vector g
% Generating from barycenter interpolation
tic
Hi=interop_bary(X,Y,xi,yi);
Hi=representant(Hi);
toc
gi=ts;
tic
Hpi=interop_bary(X,Y,xq,yq);
toc
[mm,nn]=size(Hpi);
[ii,jj,~]=find(abs(1-Hpi)<1e-5);
vv=ones(size(ii));
HpR=sparse(ii,jj,vv,mm,nn);
Hp=HpR;
ixh=indep(Hp');
ind=setdiff(1:size(Hp,1),ixh);
Hp(ind,:)=[];
gp=te*ones(size(Hp,1),1);
p.H=[Hi;Hp];
p.g=[gi;gp];
%{
% Generating from masks
jj=find(p.per1_mask);
iis=(1:length(jj))'; jjs=jj; vvs=ones(size(iis));
Hi=sparse(iis,jjs,vvs,length(iis),m*n);
gi=ts*ones(size(Hi,1),1);
jj=find(p.per12_mask);
iis=(1:length(jj))'; jjs=jj; vvs=ones(size(iis));
Hp=sparse(iis,jjs,vvs,length(iis),m*n);
gp=te*ones(size(Hp,1),1);
p.H=[Hi;Hp];
p.g=[gi;gp];
%}
%% Array u
gs=[m,n];
h=[p.dx,p.dy];
a=1.4;
bv=p.per2_time*2;
relres=1e-11;
maxit=1000;
tic
u=pdirichlet_constr(gs,h,bv,p.H,p.g,a,relres,maxit);
toc
%% Boundary conditions
p.bc=u;
%% Fuel type
p.nfuelcat=fuels2fuel(ignS.nfuel_cat);
%% Structures for computing ROS
p.ignS=ignS;
p.perS=perS;
p.perlS=perlS;
p.dynS=dynS;
%% Rate of spread
% from initial approximation
%p.R=ros_file(u,ignS,ignS,p);
%p.R(~p.vmask)=0;
% from computation using TIGN_G some time later the 2nd perimeter
%p.R=ros_file(perlS.tign_g,ignS,ignS,p);
%p.R(~p.vmask)=0;
% from ROS of the simulation at the 2nd perimeter
p.R=perS.ros;
p.R(~p.vmask)=0;
% from initial approximation and using dynamic variables
%ds=dyninterp(u,p);
%ros=ros_file(u,ignS,ds,p);
%p.ds=ds;
%p.R=ros;
%p.R(~p.vmask)=0;
%% Type of objective function and derivatives
syms x y
f=1-x*y;
p.f=f;
p.ofunc=matlabFunction(f,'Vars',[x y]);
dfdG = matlabFunction(diff(f,x),'Vars',[x y]);
dfdR = matlabFunction(diff(f,y),'Vars',[x y]);
p.dfdG=dfdG;
p.dfdR=dfdR;
%% q-norm of J
p.q=50;
%% Others
% grad_J.m
p.h=0.0001; 
% grad_min_J.m
p.stepsize=0.0001;
p.numsteps=5000;
p.max_iter=100;
p.select=@s_eno;
% linesearch.m 
p.max_step=4.;
p.nmesh=5;
p.max_depth=20;
p.min_depth=2;
p.umax=ones(m,n)*p.per2_time;
p.umin=ones(m,n)*p.per1_time;
p.bi=1:m;
p.bj=1:n;
%tign=perS.tign_g;
%tign(u<0)=nan; tign(u>p.per2_time+100)=nan;
%figure, mesh(tign'), view([0 1]);
end