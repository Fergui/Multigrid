function pdirichlet_constr_test

s = 2;     % scale up factor
n=[100*s,100*s]  % grid size
h=[1,1]
[X,Y]=ndgrid(h(1)*[0:n(1)-1],h(2)*[0:n(2)-1]);
sz = h.*n;
c = 0.5*(sz-h/2);
rmax = min(c);
cx = [c(1), c(1), c(1)]
cy = [c(2), c(2), c(2)]
r = [16, 48, 0]  
k = [100*s,200*s, 1]        % point in each perimeter last 
v = [10,40,0]                % values at each perimeter
bv= 50                    % values at the boundary
a = 1.4    % larger than values of a are needed for smooth solution
use_representant = 1

H=[];
g=[];
clf
fprintf('mesh size %i %i\n',size(X))
figure(1);
for i =1 :length(k)
    fprintf('perimeter %i points %i at %i\n',i,k(i),v(i))
    [xq,yq]=circle_points(cx(i),cy(i),r(i),k(i));  % create circles
    vq = v(i)*ones(k(i),1);
    plot3(xq,yq,vq,'.'); 
    drawnow; grid on; hold on
    Hq = interop_bary(X,Y,xq,yq);
    H = [H; Hq];
    g = [g; vq];
end
fprintf('%i variables %i constraints %i nonzeros\n',size(H'),nnz(H))

if use_representant,
    Hn=representant(H);
    [H,rows]=condense(Hn);
    g=g(rows);
end

u=pdirichlet_constr(n,h,bv,H,g,a,1e-11,1000);
figure(1);mesh(X,Y,u); xlabel('x'); ylabel('y'); zlabel('Fire arrival time'); title(['Initial approximation with alpha=',num2str(a)]);
drawnow
end

function [xq,yq]=circle_points(cx,cy,r,np)
    % create np points on circle with center c and radius r
    phi = [1:np]*2*pi/np;
    xq = cx+r*cos(phi);
    yq = cy+r*sin(phi);
end
