function [u] = primalc(S,Sinv,H,Z,f,g,varargin)
% u=primalc(S,Sinv,H,Z,f,g)
%
% Solve saddle point problem
%
%   S*u + H'*lambda = f
%   H*u             = g
% 
% with possibly singular symmetric S 
% in: 
%   Sinv    function handle, u=Sinv(r) is a solution of S*u=r if one exists
%   H       constraint matrix
%   Z       columns generate nullspace of S
%   f       column vector size number of variables
%   g       column vector size number of constraints
%   relres  max relative residual
%   maxit   max iterations
% out:
%   u       the solution

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

[m,n]=size(H);
zn = zeros(n,1);
zm = zeros(m,1);
% setup projection
% QR factorization of H'
tH=H';
R=qr(tH,0);
tR=R';
    function vn = project(u,gg)
        % orthogonal projection of u on H*u = gg
        un=H*u-gg;
        vn=u-tH*(R\(tR\un));
    end
%   PZ = @(u) u - Z*((Z'*Z)\(Z'*u));
    function v=PZ(u)
        % orthogonal projection of u on Z'*u=0
        if isempty(Z),
            v=u;
        else
            v=u-Z*((Z'*Z)\(Z'*u));
        end
    end
    % orthogonal projection on H*u = 0
    P = @(u) project(u,zm);
    
u0 = project(zn,g);    % base solution u0 that satisfies H*u0 = g
base_err=norm(H*u0-g)  
% substitution u = u0 + v
%   S*(u0 + v) + H'*lambda = f
%   H*(u0 + v)             = g
% write as
%   S*v + H'*lambda = f - S*u0
%   H*v             = 0
% which is equivalent to
%   P*S*P*v + a*(I-P)*u = P*(f - S*u0)
% where a>0 and P is the orthogonal projection on ker(H).


% set up dense problem to test vv=v
% AA=zeros(n+m);
% AA(1:n,n+1:n+m)=H';
% AA(n+1:n+m,1:n)=H;
% ei=zeros(n,1); 
% for i=1:n, ei(i)=1; AA(1:n,i)=S(ei); ei(i)=0; end
% bb =[f - AA(1:n,1:n)*u0; zeros(m,1)];
% vl=AA\bb;
% vv=vl(1:n);

r = P(f-S(u0));
reg=1;
A = @(u) P(S(P(u)))+reg*(u-P(u));
M = @(r) P(PZ(Sinv(PZ(P(r)))));
[v,FLAG,RELRES,ITER,resvec]=pcg(A,r,relres,maxit,M);
ITER,RELRES
u=u0+P(v);
res1 = norm(P(S(u) - f));
res2 = norm(H*u - g);

fprintf('primalc: residual %g %g\n',res1,res2)

end
