function plot_sol(u,H,g)
%plot_sol(u,H,g)
% Plot the fire arrival time and the constraints in a 3D contour
%in
%   u   fire arrival time over the domain
%   H   matrix of the constraints Hu=g
%   g   right hand side of the constraints Hu=g

col=['-g','-b','-r'];
uu=unique(g);
for k=1:length(uu)
    gs=g==uu(k);
    us=H'*gs;
    us=reshape(us,size(u));
    contour3(flip(uu(k)*(us>0)),[uu(k),uu(k)],col(k)), hold on
end

contour3(flip(u),[uu(1):(uu(end)-uu(1))/10:uu(end)],'k')
%contour3(flip(u),200,'k')

end

