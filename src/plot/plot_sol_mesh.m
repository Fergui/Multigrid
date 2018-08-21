function plot_sol_mesh(X,Y,u,H,g)
%plot_sol_mesh(u,H,g)
% Plot the fire arrival time and the constraints in a 3D contour and mesh
%in
%   X       matrix of X coordinates (longitude), or number of points in X direction
%   Y       matrix of Y coordinates (latitude), or number of points in Y direction
%   u   fire arrival time over the domain
%   H   matrix of the constraints Hu=g
%   g   right hand side of the constraints Hu=g

[nx,ny]=size(X);
if nx==1,
    nx=X;
    ny=Y;
    [X,Y]=ndgrid(1:nx,1:ny);
end

mesh(X,Y,u), hold on

col=[255,0,0;255,165,0;0,100,255]/255;
uu=unique(g);
for k=1:length(uu)
    gs=g==uu(k);
    us=H'*gs;
    us=reshape(us,size(u));
    %contour3(flip((uu(2)*2)*(us>0)),[uu(2)*2,uu(2)*2],col(k)), hold on
    contour3(X,Y,(uu(2))*(us>0),[uu(2),uu(2)],'Color',col(k,:)), hold on
end

end