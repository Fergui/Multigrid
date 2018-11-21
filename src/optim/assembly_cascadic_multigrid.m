function [mg,levels]=assembly_cascadic_multigrid(p,times)
%% Setup
if ~isfield(p,'dynR')
    p.dynR=0;
end
if ~isfield(p,'interp')
    p.interp='cubic';
end
if ~isfield(p,'coarse')
    p.coarse=0;
end
if ~isfield(p,'max_iter')
    p.max_iter=20;
end
% determine ratio and finest grid mesh
ratio=2;
sm=findmgfit(p.m);
sn=findmgfit(p.n);
m=length(sm);
n=length(sn);
% determine number of levels
mm=m; nn=n;
for levels=1:100,
    mm=(mm-1)/ratio+1; 
    nn=(nn-1)/ratio+1;
    if mm~=round(mm) || nn~=round(nn) || mm<4 || nn<4
        break
    end
end
levels=4;

%% Initialize level 1
fprintf('assembly_cascadic_multigrid.m: initializing level 1.\n');
mg(1)=p;
% Dimensions and resolutions
mg(1).m=m;
mg(1).n=n;
mg(1).dx=p.dx;
mg(1).dy=p.dy;
% Grid and indexes
mg(1).ii=1:mg(1).m;
mg(1).jj=1:mg(1).n;
[mg(1).II,mg(1).JJ]=meshgrid(mg(1).ii,mg(1).jj);
mg(1).iic=1:ratio:mg(1).m;
mg(1).jjc=1:ratio:mg(1).n;
mg(1).X=p.X(sm,sn);
mg(1).Y=p.Y(sm,sn);
mg(1).Xc=mg(1).X(mg(1).iic,mg(1).jjc);
mg(1).Yc=mg(1).Y(mg(1).iic,mg(1).jjc);
% H and g for the linear system of constraints
[Hs,gs]=interpolation(mg(1).X,mg(1).Y,p.shapes,times);
Hc=cat(1,Hs{:});
[H,rows]=condence(Hc);
[mg(1).H,kk]=repre(H);
gc=cat(1,gs{:});
rows=rows(kk);
g=gc(rows);
mg(1).g=g;
% Other configuration parameters
mg(1).vmask=p.vmask(sm,sn);
% ROS dynamic important variables
if p.dynR
    ffuel=fields(p.nfuelcat);
    for k=1:length(ffuel)
        mg(1).nfuelcat.(ffuel{k})=p.nfuelcat.(ffuel{k})(sm,sn);
    end
    mg(1).ignS.dzdxf=p.ignS.dzdxf(sm,sn);
    mg(1).ignS.dzdyf=p.ignS.dzdyf(sm,sn);
    mg(1).dynS.simt=p.dynS.simt;
    mg(1).dynS.uf=p.dynS.uf(sm,sn,:);
    mg(1).dynS.vf=p.dynS.vf(sm,sn,:);
    mg(1).dynS.fmc_g=p.dynS.fmc_g(sm,sn,:);
else
    R=p.R(sm,sn);
    mg(1).R=R;
end

%% Prepare all the levels
for i=2:levels
    fprintf('assembly_cascadic_multigrid.m: prepare level %d.\n',i);
    mg(i)=mg(1);
    % Dimensions and resolutions
    mg(i).m=length(mg(i-1).iic);
    mg(i).n=length(mg(i-1).jjc);
    mg(i).dx=mg(i-1).dx*ratio;
    mg(i).dy=mg(i-1).dy*ratio;
    % Grid and indexes
    mg(i).ii=mg(i-1).iic; 
    mg(i).jj=mg(i-1).jjc;
    [mg(i).II,mg(i).JJ]=meshgrid(1:2^(i-1):mg(1).m,1:2^(i-1):mg(1).n);
    mg(i).iic=1:ratio:mg(i).m;
    mg(i).jjc=1:ratio:mg(i).n;
    mg(i).X=mg(i-1).Xc;
    mg(i).Y=mg(i-1).Yc;
    mg(i).Xc=mg(i).X(mg(i).iic,mg(i).jjc);
    mg(i).Yc=mg(i).Y(mg(i).iic,mg(i).jjc);
    % H and g for the linear system of constraints
    [Hs,gs]=interpolation(mg(i).X,mg(i).Y,p.shapes,times);
    Hc=cat(1,Hs{:});
    [H,rows]=condence(Hc);
    gc=cat(1,gs{:});
    [mg(i).H,kk]=repre(H);
    rows=rows(kk);
    mg(i).g=gc(rows);
    % Rate of spread
    if mg(i).coarse
        mg(i).R=coarsening(mg(i-1).R,2);
        mg(i).R(~mg(i).vmask)=0;
    else
        mg(i).R=mg(i-1).R(mg(i).ii,mg(i).jj);
    end
    % Other configuration parameters
    mg(i).vmask=mg(i-1).vmask(mg(i).ii,mg(i).jj);
    %mg(i).vmask=(coarsening(mg(i-1).vmask,2)>0);
    % ROS dynamic important variables
    if p.dynR
        ffuel=fields(p.nfuelcat);
        for k=1:length(ffuel)
            mg(i).nfuelcat.(ffuel{k})=mg(i-1).nfuelcat.(ffuel{k})(mg(i).ii,mg(i).jj);
        end
        mg(i).ignS.dzdxf=mg(i-1).ignS.dzdxf(mg(i).ii,mg(i).jj);
        mg(i).ignS.dzdyf=mg(i-1).ignS.dzdyf(mg(i).ii,mg(i).jj);
        mg(i).dynS.simt=mg(i-1).dynS.simt;
        mg(i).dynS.uf=mg(i-1).dynS.uf(mg(i).ii,mg(i).jj,:);
        mg(i).dynS.vf=mg(i-1).dynS.vf(mg(i).ii,mg(i).jj,:);
        mg(i).dynS.fmc_g=mg(i-1).dynS.fmc_g(mg(i).ii,mg(i).jj,:);
    end
end

%% Last level initialization
fprintf('assembly_cascadic_multigrid.m: last level initialization.\n');
% Initial fire arrival time using Dirichlet boundary conditions
uu=unique(mg(end).g);
gs=[mg(end).m,mg(end).n];
h=[mg(end).dx,mg(end).dy];
a=1;
bv=uu(end)*1.05;
relres=1e-11;
maxit=1000;
tic
u=pdirichlet_constr(gs,h,bv,mg(end).H,mg(end).g,a,relres,maxit);
toc
u=max(0,u);
mg(end).u=u;
% Boundary conditions
mg(end).bc=u;

end
