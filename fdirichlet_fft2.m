function U = fdirichlet_fft2(r,h,f)
% v=fdirichlet_fft2(r,h,f)
% Multiply by a function of minus Laplace operator with dirichlet zero
% boundary conditions by sine FFT
%
% in:
%   r   2d gridded array
%   h   vector size 2: grid spacing
%   f   function, must work entry by entry for matrix matrix argument 
% out:
%   u = f(-d^2/dx^2 - d^2/dy^2)*r

% with [N1,N2]=size(r) on the square [0 (N1+1)*h(1)] times [0 (N2+1)*h(2)]
% return values at nodes [i1*h(1),i2*h(2)], i1=1:N1, i2=1:N2
% with one layer reflected about the boundary

[N1,N2] = size(r);
% eigenvalue terms
e1=poisson_1d_eig(N1,h(1));
e2=poisson_1d_eig(N2,h(2));
e=zeros(N1,N2);
for i2=1:N2
    for i1=1:N1
        e(i1,i2)=e1(i1)+e2(i2);
    end
end
fe=f(e);
% VV=(2/(n+1)*ones(n,n)./(V'(ones(1,n)+ones(1,n)*V');
U=dst2(r);
U=U.*fe;
U=dst2(U)*4/((N1+1)*(N2+1)); % scale for nonunitary DST2



