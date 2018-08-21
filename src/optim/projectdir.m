function d = projectdir(d,H)
%v = projectgrad(u,H)
% Computes the projected direction for projected method
%in
%   dir  search direction
%   H    sparse matrix, interpolation operator (constraint Hd=0)
%out
%   d    direction to search for in the linesearch

[m,n]=size(d);
% QR factorization of H'
tH=H';
R=qr(tH);
tR=R';
% Projecting the array u
dn=H*d(:);
if nnz(dn)
    d=d(:)-tH*(R\(tR\dn));
    d=reshape(d,m,n);
end

end