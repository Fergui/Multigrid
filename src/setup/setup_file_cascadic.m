function p = setup_file_cascadic(file,perF,dynR)
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
if dynR
    dynS=nc2struct(file,{'UF','VF','FMC_G'},{});
    dynS.simt=simt;
end
ts=simt(p.ignF);
te=simt(p.perF);
%% Dimensions and spacings
[m,n]=size(ignS.fire_area);
sr_x=round(size(ignS.fxlong,1)/size(ignS.xlong,1),1);
sr_y=round(size(ignS.fxlat,1)/size(ignS.xlat,1),1);
p.m=m;
p.n=n;
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
%% Fuel type
p.nfuelcat=fuels2fuel(ignS.nfuel_cat);
%% Structures for computing ROS
p.ignS=ignS;
p.perS=perS;
%% Rate of spread
if dynR
    p.dynS=dynS;
    ds=dyninterp(p.perS.tign_g,p);
    ros=ros_file(p.perS.tign_g,ignS,ds,p);
    p.R=ros;
    p.R(~p.vmask)=0;
else
    p.R=p.perS.ros;
    p.R(~p.vmask)=0;
end
%% Type of objective function and derivatives
syms x y
f=1-x*y;
p.f=f;
p.ofunc=matlabFunction(f,'Vars',[x y]);
%% q-norm of J
p.q=4;
%% Others
p.select=@s_eno;
p.max_step=1.0;
p.nmesh=5;
p.max_depth=2;
p.min_depth=1;
p.exp='file';
p.penalty=1;
p.X=X;
p.Y=Y;
p.dynR=dynR;
end
