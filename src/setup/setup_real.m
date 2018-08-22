function [u,s,p] = setup_real(matfile)
% Call:
% [u,s,p] = setup_real(matfile)
%
% Example: 
% [u,p,s] = setup_real('in/lasconchas.mat');
%
% Inputs: 
%   matfile   Matlab file with an structure p on it containing all the 
%             necessary information from wrfout and kml 
% Outputs:
%   u        struct of level sets initialized using an approximation for
%            each perimeters case
%   s        structure with:
%               dx, dy      fire mesh spacing
%               H           structure of interpolation operator matrices (one for each perimeter)
%               g           structure of right sides of Hu=b for each H
%               M           structure of masks for each perimeter
%               bc          structure of boundary conditions for each case = fixed values of the level set function where mask is false
%               nfuelcat    structure with all the fuel data necessary to compute ROS
%               R           matrix, rate of spread, on same nodes as u
%               f           string with the objective function formula
%               ofunc       matlab function, objective function comparing x=||grad u||^2 and y=R^2 such that xy=1
%               dfdG        matlab function, partial derivative of ofunc to respect to x=||grad u||^2
%               dfdG        matlab function, partial derivative of ofunc to respect to y=R^2
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
%				bi			indeces to compute the first objective function (coordinate x)
%				bj			indeces to compute the first objective function (coordinate y)
%               exp         experiment type, string 'real'
%   p       structure with:
%               sdates      simulation dates
%               stimes      simulation times from the simulation start
%               tig         ignition time from the simulation start
%               kml         structure with the real perimeters information
%               pdates      perimeter dates
%               ptimes      perimeter times from the simulation start
%               sframes     simulation frames where the perimeters are from
%               ignS        structure with all the important variables from the simulation in the ignition time
%               perS        structure with all the important variables from the simulation in the perimeter times
%               perlS       structure with all the important variables from the simulation a posteriori of the perimeter times
%               dynS        structure with the dynamic variables components winds (uf,vf) and fuel moisture (fmc_g) at all the time steps
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Take data from files
load(matfile);

fprintf('Computing necessary mesh variables...\n');
%% Fire mesh refinement
sr_x=round(size(p.ignS.fxlong,1)/size(p.ignS.xlong,1),1);
sr_y=round(size(p.ignS.fxlat,1)/size(p.ignS.xlat,1),1);

%% Number of points, spacing and mesh
[m,n]=size(p.ignS.fxlat);
X=p.ignS.fxlong;
Y=p.ignS.fxlat;
s.dx=p.ignS.dx/sr_x;
s.dy=p.ignS.dy/sr_y;
np=length(p.perS);
minp=500; % Minim of points in the perimeters
times=[p.tig,p.ptimes(1:np)];

%% Defining the shape points
xq=cell(np+1,1);
yq=cell(np+1,1);
% Ignition point
xq{1}=p.ignS.coordinates(1);
yq{1}=p.ignS.coordinates(2);
% Computing perimeters
for k=1:np
    xlon=p.kml(k).Lon;
    xlat=p.kml(k).Lat;
    for kk=1:length(xlon)
        lon=xlon{kk};
        lat=xlat{kk};
        if length(lon)>minp        
            xx=lon(1:end-2);
            yy=lat(1:end-2);
            [lons,lats]=addshape(X,Y,xx,yy);
            xq{k+1}=[xq{k+1};lons];
            yq{k+1}=[yq{k+1};lats];
        end
    end
end

fprintf('Computing interpolation operator H and right side g of Hu=g...\n');
%% Interpolate operator H and g such that Hu=g
H=cell(np+1,1);
g=cell(np+1,1);
for k=1:np+1
    fprintf('Interpolation of perimeter: %d/%d\n',k,np+1);
    tic
    I=interop_bary(X,Y,xq{k},yq{k});
    toc
    tic
    H{k}=condense(I);
    toc
    g{k}=times(k)*ones(size(H{k},1),1);
end
s.H=H;
s.g=g;

fprintf('Computing approximation of the fire arrival time...\n');
%% Masks
M=cell(np,1);
v=cell(np,1);
for k=1:np
    fprintf('Mask between ignition point and consecutive perimeters: %d/%d\n',k,np);
    Hm=[H{1};H{k+1}];
    gm=[g{1};g{k+1}];
    Rm=representant(Hm);
    [Cm,rm]=condense(Rm);
    dm=gm(rm);
    tic
    [M{k},v{k}]=per2mask(Cm,dm,m,n,s.dx,s.dy);
    toc
end
s.M=M;

%% Structure of arrays u
u=cell(np,1);
for k=1:np
    fprintf('Condense generating artificial perimeter: %d/%d...\n',k,np);
    Hm=[H{k};H{k+1}]; 
    gm=[g{k};g{k+1}];
    tic
    [Hc,gc]=condense_arti(Hm,gm,v{k});
    toc
    uu=unique(gc);
    gs=[m,n];
    h=[s.dx,s.dy];
    a=1.4;
    bv=uu(end)*2;
    relres=1e-11;
    maxit=1000;
    tic
    u{k}=pdirichlet_constr(gs,h,bv,Hc,gc,a,relres,maxit);
    toc
    u{k}=max(0,u{k});
end

fprintf('Computing fuel variables...\n'); 
%% Fuel type
tic
s.nfuelcat=fuels2fuel(p.ignS.nfuel_cat);
toc

fprintf('Computing ROS...\n'); 
%% Rate of spread
% from fire arrival time initial approximation
R=cell(np,1);
R{1}=ros_file(u{1},p.ignS,p.ignS,s);
for k=2:np
    R{k}=ros_file(u{k},p.ignS,p.perS(k-1),s);
end
s.Ri=R;
% from fire arrival time at the end of the simulation and winds from
% ignition frame of wrfout
T=cell(np,1);
for k=1:np
    T{k}=p.perlS(k).tign_g;
end
s.T=T;
R=cell(np,1);
R{1}=ros_file(T{1},p.ignS,p.ignS,s);
for k=2:np
    R{k}=ros_file(T{k},p.ignS,p.perS(k-1),s);
end
s.Rd=R;
% from ROS at the perimeter time
R=cell(np,1);
for k=1:np
   R{k}=p.perS(k).ros; 
end
s.R=R;

fprintf('Defining final variables...\n');
%% Boundary conditions
s.bc=u;
%% Type of objective function and derivatives
syms x y
f=x*y-1;
s.f=f;
s.ofunc=matlabFunction(f,'Vars',[x y]);
dfdG = matlabFunction(diff(f,x),'Vars',[x y]);
dfdR = matlabFunction(diff(f,y),'Vars',[x y]);
s.dfdG=dfdG;
s.dfdR=dfdR;
%% q-norm of J
s.q=4;
%% Others
% grad_J.m
s.h=0.0001; 
% grad_min_J.m
s.stepsize=0.0001;
s.numsteps=5000;
s.max_iter=100;
s.select=@s_eno;
% linesearch.m 
s.max_step=4.;
s.nmesh=5;
s.max_depth=20;
s.min_depth=2;
umin=cell(np,1); % Min and max
umax=cell(np,1);
for k=1:np
    umin{k}=ones(m,n)*times(k);
    umax{k}=ones(m,n)*times(k+1);
end
s.umin=umin; 
s.umax=umax;
s.bi=1:m; % First patch generation
s.bj=1:n;
end
