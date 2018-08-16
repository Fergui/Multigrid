% Example

% Taking the shape file
s=kml2struct('data/doc.kml');

% Using grid from wrfout
m=2580;n=2580;
ulon=-106.0347;ulat=36.1729;llon=-106.7677;llat=35.5790;

% Defining important data computing H
x=linspace(llon,ulon,m);
y=linspace(llat,ulat,n);
[X,Y]=meshgrid(x,y);
xq=s(1).Lon{3};
xq=xq(1:end-2);
yq=s(1).Lat{3};
yq=yq(1:end-2);

tic
H=interop_bary(X,Y,xq,yq);
toc

g=909*ones(1,size(H,1));

figure, plot_constr_scatter(X,Y,H,g);
