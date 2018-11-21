function p=setup_circles_cascadic(m,n,dx,dy)
%% Dimensions and spacings
p.m=m;
p.n=n;
p.dx=dx;
p.dy=dy;
%% Mask perimeters
x0=dx*(m-1)/2;
y0=dy*(n-1)/2;
d1=14;
d2=42;
per1=false(m,n);
per2=false(m,n);
for i=1:m
    for j=1:n
        x=dx*(i-1);
        y=dy*(j-1);
        r=(x-x0)^2+(y-y0)^2;
        if r<=d1^2
            per1(i,j)=1;
        end
        if r>=d2^2
            per2(i,j)=1;
        end
    end
end
p.per1_mask=per1;
p.per2_mask=per2;
p.per12_mask=per2;
p.per0_time=0;
p.per1_time=19;
p=setup_masks(p);
%% Rate of spread
p.R=zeros(m,n);
T1=p.per1_time;
R1=1;
R2=2;
D1=d1;
D2=d1;
T2=T1+D1/R1+D2/R2;
p.per2_time=T2;
for i=1:m
    for j=1:n
        x=dx*(i-1);
        y=dy*(j-1);
        r=(x-x0)^2+(y-y0)^2;
        if p.mask(i,j)
            if r<=(d1+D1)^2
                p.R(i,j)=R1;
            end
            if r>(d2-D2)^2
                p.R(i,j)=R2;
            end
        end
    end
end
p.R(p.per1_mask)=14/19;
p.mask(p.per1_mask)=1;
p.vmask=p.mask;
%% Shape points
[xq1,yq1]=circle_points(x0,y0,d1,150);
[xq2,yq2]=circle_points(x0,y0,d2,250);
[X,Y]=meshgrid(0:m-1,0:n-1);
[xq1n,yq1n]=addshape(X,Y,xq1',yq1');
[xq2n,yq2n]=addshape(X,Y,xq2',yq2');
p.shapes(1).x=x0;
p.shapes(1).y=y0;
p.shapes(2).x=xq1n;
p.shapes(2).y=yq1n;
p.shapes(3).x=xq2n;
p.shapes(3).y=yq2n;
%% Type of objective function
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
p.max_depth=10;
p.min_depth=2;
p.exp='ideal';
p.penalty=1;
p.X=X;
p.Y=Y;
p.coarse=0;
end

function [xq,yq]=circle_points(cx,cy,r,np)
    % create np points on circle with center c and radius r
    phi = [1:np]*2*pi/np;
    xq = cx+r*cos(phi);
    yq = cy+r*sin(phi);
end