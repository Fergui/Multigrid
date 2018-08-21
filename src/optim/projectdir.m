function d = projectdir(d,H)
% Call:
% v = projectgrad(u,H)
%
% Description:
% Projects the search direction into Hd=0
%
% Inputs
%   d    search direction
%   H    sparse matrix, interpolation operator (constraint Hd=0)
% Outputs:
%   d    d search direction projected to Hd=0
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

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
