function [u,p] = setup_ideal(m,n,dx,dy,upwind,type)
% Call:
% [u,p] = setup_ideal(m,n,dx,dy,upwind)
%
% Examples:
%[u,p] = setup_ideal(50,50,1,1,@s_eno,'s');
%[u,p] = setup_ideal(50,50,1,1,@s_eno,'c');
%
% Inputs:
%   m       x dimension
%   n       y dimension
%   dx      x resolution
%   dy      y resolution
%   upwind  upwind method (example: @s_eno) 
%   type    Type of experiment, possible values:
%         s - Squares experiment.
%         c - Circles experiment.
% Outputs:
%   u     level set initialized using an approximation
%   p        p structure with:
%               dx, dy      mesh spacing
%               per1_mask   matrix, mask of the first perimeter
%               per12_mask  matrix, mask of the second perimeter line
%               per2_mask   matrix, mask of the second perimeter
%               per1_time   fire arrival time of the first perimeter
%               per2_time   fire arrival time of the second perimeter
%               mask        matrix, true where the values of level set function have to change
%               vmask       matrix, true where the values of level set function can change
%               bc          boundary conditions = fixed values of the level set function where mask is false
%               R           matrix, rate of spread, on same nodes as u
%               H           sparse matrix, interpolation operator for all the perimeters as a constraints
%               g           vector, right side of the constraints (Hu=g)
%               bc          boundary conditions for each case = fixed values of the level set function where mask is false
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
%               bi			indeces to compute the first objective function (coordinate x)
%               bj			indeces to compute the first objective function (coordinate y)
%               exp         experiment type, string 'ideal'
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Perimeters
% Squares
if type=='s'
    per2=false(m,n);
    per2(1,:)=1;per2(end,:)=1;per2(:,1)=1;per2(:,end)=1;
    per2(2,:)=1;per2(end-1,:)=1;per2(:,2)=1;per2(:,end-1)=1;
    per2(3,3)=1;per2(end-2,3)=1;per2(3,end-2)=1;per2(end-2,end-2)=1;
    per1=false(m,n);
    per1(23:27,23:27)=1;
% Circles 
elseif type=='c'
    x0=dx*(m-1)/2;
    y0=dy*(n-1)/2;
    d1=16;
    d2=48;
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
else
    error('Bad type of experiment used');
end
%% Spacing
p.dx=dx;
p.dy=dy;
%% Masks
p.per1_mask=per1;
p.per2_mask=per2;
p.per12_mask=per2;
p.per1_time=16;
% Squares
if type=='s'
    p.per2_time=10;
end
p=setup_masks(p);
%% Rate of spread (R)
% Squares
if type=='s'
    p.R=ones(m,n)*5;
    p.R(:,7:12)=7;
    p.R(:,13:15)=10;
% Circles
elseif type=='c'
    p.R=zeros(m,n);
    T1=p.per1_time;
    R1=1;
    R2=2;
    D1=16;
    D2=16;
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
    %% Creating shape
    [xq1,yq1]=circle_points(x0,y0,d1,100);
    [xq2,yq2]=circle_points(x0,y0,d2,200);
    x=0:m-1;
    y=0:n-1;
    [X,Y]=meshgrid(x,y);
    [xq1n,yq1n]=addshape(X,Y,xq1',yq1');
    [xq2n,yq2n]=addshape(X,Y,xq2',yq2');
    H1=interop_bary(X,Y,xq1n,yq1n);
    H2=interop_bary(X,Y,xq2n,yq2n);
    H=[H1;H2];
    [p.H,rows]=condense(H);
    b1=p.per1_time*ones(size(H1,1),1);
    b2=p.per2_time*ones(size(H2,1),1);
    g=[b1;b2];
    p.g=g(rows);
end
%% Initialization of u using Dirichlet boundary conditions
uu=unique(p.g);
gs=[m,n];
h=[p.dx,p.dy];
a=1.4;
bv=uu(end);
relres=1e-11;
maxit=1000;
tic
u=pdirichlet_constr(gs,h,bv,p.H,p.g,a,relres,maxit);
toc
u=max(0,u);
%% Boundary conditions
p.bc=u;
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
p.q=10;
%% Others
% grad_J.m
p.h=0.0001; 
% grad_min_J.m
p.stepsize=0.0001;
p.numsteps=5000;
p.max_iter=25;
p.select=upwind;
% linesearch.m 
p.max_step=4.0;
p.nmesh=5;
p.max_depth=20;
p.min_depth=2;
p.umax=ones(m,n)*p.per2_time;
p.umin=ones(m,n)*p.per1_time;
p.bi=1:m;
p.bj=1:n;
p.exp='ideal';
end

function [xq,yq]=circle_points(cx,cy,r,np)
    % create np points on circle with center c and radius r
    phi = [1:np]*2*pi/np;
    xq = cx+r*cos(phi);
    yq = cy+r*sin(phi);
end
