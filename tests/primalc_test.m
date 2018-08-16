function primalc_test 
m=300; n=3;
S = randn(m); S=S*S';
Z = [];
Sinv = pinv(S);
f = randn(m,1);
g = randn(n,1);
H = sparse(randn(n,m));
[u] = primalc(@(r) S*r, @(x) Sinv*x ,H,Z,f,g);
end