function [u,p] = setup_file(file,perF)
% Call:
% [u,p] = setup_file(file,perT)
%
% Example: 
% [u,p] = setup_file('wrfout',50);
%
% Description:
% Set up fire arrival time u and matlab struct p for the fire arrival time interpolation method between 
% the ignition point and a perimeter at frame perF in the case of WRF-SFIRE ideal simulations.
%
% Inputs: 
%   file     NetCDF WRF-SFIRE output file
%   perF     wrfout second perimeter frame
% Outputs:
%   u        level set initialized using an approximation
%   p        Matlab structure with:
%               perF        wrfout perimeter frame
%               lastF       wrfout last frame
%               ignF        wrfout ignition frame
%               dx,dy       Fire mesh spacing
%               per1_mask   Matrix, mask of the first perimeter
%               per12_mask  Matrix, mask of the second perimeter line
%               per2_mask   Matrix, mask of the second perimeter
%               per1_time   Fire arrival time of the first perimeter
%               per2_time   Fire arrival time of the second perimeter
%               mask        Matrix, true where the values of level set function have to change
%               vmask       Matrix, true where the values of level set function can change
%               H           Interpolation operator matrix of the ignition point and perimeter at frame perF
%               g           Right hand side of Hu=g, i.e. fire arrival times at the constraints
%               bc          boundary conditions = fixed values of the level set function where mask is false
%               nfuelcat    Structure with all the fuel data necessary to compute ROS
%               ignS        Matlab structure with:
%                       timestep        Number of frame at ignition time
%                       filename        wrfout file name
%                       fxlong,fxlat    Fire longitude and latitude coordinates
%                       xlong,xlat      Atmospheric longitude and latitude coordinates
%                       fire_area       Mask of the fire area at the ignition time
%                       nfuel_cat       Fuel types all over the domain (13 Rothermel's classification)
%                       dzdxf,dzdyf     Slope components at ignition time
%                       uf,vf           Wind components at ignition time
%                       fmc_g           Fuel moisture content at ignition time
%                       dx,dy           Resolutions of the atmospheric grid
%               perS        Matlab structure with:
%                       timestep        Number of frame at perimeter time = p.perF
%                       filename        wrfout file name
%                       fire_area       Mask of the fire area at the perimeter time
%                       tign_g          Fire arrival time at the perimeter time
%                       ros             Rate of spread at the perimeter time
%               dynS        Matlab structure with:
%                       timestep        =0, because the data is at all the timesteps
%                       filename        wrfout file name
%                       uf,vf           Wind components at all the timesteps
%                       fmc_g           Fuel moisture at all the timesteps
%                       simt            Seconds from the start of the simulation of all the timesteps
%               R           Matrix, rate of spread, on same nodes as u
%               f           
%               ofunc       Matlab function, objective function comparing x=||grad u||^2 and y=R^2 such that xy=1
%               dfdG        Matlab function, partial derivative of ofunc to respect to x=||grad u||^2
%               dfdG        Matlab function, partial derivative of ofunc to respect to y=R^2
%               q           q norm of the computation of J
%               h           The stepsize to compute the gradient
%               stepsize    Step size for minimization
%               numsteps    Number of steps to try
%               max_iter    Maximum number of iterations
%               select      Handle to upwinding function
%               max_step    Max size of the steps to search
%               nmesh       Number of mesh points in each search
%               max_depth   Max number of searchs
%               min_depth   Min number of searchs
%               umax        Array, maximal value of u
%               umin        Array, minimal value of u
%				bi			Indeces to compute the first objective function (coordinate x)
%				bj			Indeces to compute the first objective function (coordinate y)
%               exp         experiment type, string 'file'
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

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
dynS=nc2struct(file,{'UF','VF','FMC_G'},{});
dynS.simt=simt;
ts=simt(p.ignF);
te=simt(p.perF);
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
%% Masks
p.per1_mask=logical(per1);
p.per12_mask=per12;
p.per2_mask=per2;
p.per1_time=ts;
p.per2_time=te;
p=setup_masks(p);
%% Shape points
X=ignS.fxlong;
Y=ignS.fxlat;
xi=xx*p.dx;
yi=yy*p.dy;
p.shapes(1).x=xi;
p.shapes(1).y=yi;
xq=X(per12);
yq=Y(per12);
p.shapes(2).x=xq;
p.shapes(2).y=yq;
%% Matrix H and vector g
[Hs,gs]=interpolation(X,Y,p.shapes,[p.per1_time,p.per2_time]);
H=[Hs{1};Hs{2}];
[p.H,rows]=condence(H);
g=[gs{1};gs{2}];
p.g=g(rows);
%{
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
%}
%% Array u
gs=[m,n];
h=[p.dx,p.dy];
a=1.4;
bv=p.per2_time*5;
relres=1e-11;
maxit=1000;
tic
u=pdirichlet_constr(gs,h,bv,p.H,p.g,a,relres,maxit);
toc
u=max(0,u);
%% Boundary conditions
p.bc=u;
%% Fuel type
p.nfuelcat=fuels2fuel(ignS.nfuel_cat);
%% Structures for computing ROS
p.ignS=ignS;
p.perS=perS;
p.dynS=dynS;
%% Rate of spread
ds=dyninterp(u,p);
ros=ros_file(u,ignS,ds,p);
p.R=ros;
p.R(~p.vmask)=0;
%p.R=p.perS.ros;
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
p.q=4;
%% Others
% grad_J.m
p.h=0.0001; 
% grad_min_J.m
p.stepsize=0.0001;
p.numsteps=5000;
p.max_iter=100;
p.select=@s_eno;
% linesearch.m 
p.max_step=1.0;
p.nmesh=5;
p.max_depth=2;
p.min_depth=1;
p.mcycle=4;
p.umax=ones(m,n)*p.per2_time;
p.umin=ones(m,n)*p.per1_time;
p.bi=1:m;
p.bj=1:n;
p.exp='file';
p.penalty=1;
p.X=X;
p.Y=Y;
end
