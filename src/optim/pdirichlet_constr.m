function u=pdirichlet_constr(n,h,bv,H,g,a,varargin)
% Call:
% u=pdirichlet_constr(n,h,bv,H,g,a)
%
% Description:
% Solve   (-L)^a * u = 0 on rectangle 
%         u          = bv on the boundary, subject to   
%         H*u        = g
% where L = d^2/dx^2 + d^2/dy^2 is the Laplace operator
%
% Inputs:
%   n   vector size 2, the size of the rectangle in nodes in x and y
%   h   vector size 2: mesh step
%   bv   boundary value, scalar
%   H   the constraint matrix, multiply gridded array as a single column
%   g   column, contraint right hand sides 
%   a   exponent, should be >1 for best results
%
% Outputs:
%   u   the solution, matrix size n(1) by n(2)
%
% Jan Mandel, 2018
%-------------------------------------------------------------------------

if ~isempty(varargin)
    relres=varargin{1};
else
    relres=1e-6;
end
if length(varargin)>1
    maxit=varargin{2};
else
    maxit=500;
end

check_args;

% remove redundant constraints
[k0,m]=size(H);  % total number of variables and constraints
nn=prod(n);
fprintf('pdirichlet_constr: domain %g by %g variables %g exponent %g\n',n,nn,a)
ix = indep(H');
HH=H(ix,:);
g=g(ix);
[k,m]=size(HH);  % total number of variables and constraints
fprintf('pdirichlet_constr: constraints %g linearly independent %g\n',k0,k)

% substituting  u = u0 + b, b = constant function value bv 
%        (-L)^a * (u0 + b) = 0 
%         u0               = 0  on the boundary, subject to  
%         H*(u0 + b)       = g
% gives using that (-L)^a * constant function = 0
%        (-L)^a * u0  = 0 
%         u0          = 0  on the boundary, subject to  
%         H*u0        = g - H*b


Sinv = @(r) Smult(r,-a);
S = @(r) Smult(r,a);
Z = [];
gg = g - HH*bv*ones(nn,1);  % or sum(HH,2)
ff = -.1*ones(nn,1); % cannot -S(b(:)) because S has zero boundary conditions inside
uu=primalc(S,Sinv,HH,Z,ff,gg,relres,maxit);

u =reshape(uu,n) + bv;

    function v=Smult(u,e)
        % v=(-Laplace)^2*u
        % input and output as single column
        v = fdirichlet_fft2(reshape(u,n),h, @(t) t.^(e));
        v = v(:);
    end

    function check_args
        if ~ismatrix(H) 
            error('H must be matrix')
        end
        if ~issparse(H)
            warning('For best results, H should be sparse')
        end
        if prod(n) ~= size(H,2)
            error('incompatible domain size n and size of H')
        end
        if length(g) ~= size(H,1)
            error('incompatible sizes of g and H')
        end
        if ~isvector(n) || length(n) ~= 2
            error('n must be vector length 2')
        end
        if ~isvector(h) || length(h) ~= 2
            error('h must be vector length 2')
        end
        if ~isscalar(a)
            error('a must be scalar')
        end
        if ~isscalar(bv)
            error('b must be scalar')
        end
        
    end

end
