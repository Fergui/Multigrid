function ix=indep(A,tol)
% ix=indep(A)
% ix=indep(A,tol)
% select linearly independent columns in matrix
% in
%   A matrix
%   t tolerance to identify nonzeros, default eps*max(size(A))
% out
%   ix indices of a maximal set of linearly independent columns

if ~exist('tol','var')
    tol=max(size(A))*eps;
end
R = qr(A,0);
% for each row j of R, find the index of the first nonzero. 
% Then the column is linearly independent of previous columns
[m,n]=size(A);
ix=zeros(1,n);
t=tol*max(abs(R(:)));
k=0;
for j = 1:n
    if(abs(R(k+1,j)) > t)
        k=k+1;
        ix(k)=j;
    end
end
ix=ix(1:k);
end
