function plot_sol(u,H,g)
% Call:
% plot_sol(u,H,g)
%
% Description:
% Plot the fire arrival time and the constraints in a 3D contour
%
% Inputs:
%   u   fire arrival time over the domain
%   H   matrix of the constraints Hu=g
%   g   right hand side of the constraints Hu=g
% Outputs:
%   The result is a 3D plot of the solution u with the constraints Hu=g.
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

col=['-g','-b','-r'];
uu=unique(g);
for k=1:length(uu)
    gs=g==uu(k);
    us=H'*gs;
    us=reshape(us,size(u));
    contour3(flip(uu(k)*(us>0)),[uu(k),uu(k)],col(k)), hold on
end

contour3(flip(u),[uu(1):(uu(end)-uu(1))/10:uu(end)],'k')

end

