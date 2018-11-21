function H = representant(A)
% Call:
% H = representant(A)
%
% Description:
% Take one unique non-zero value (representative) for each row of A
%
% Inputs: 
%   A   Matrix to take the representants for rows
% Outputs:
%   H   Final matrix with 1's in each representant of each row
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

% Dimensions of the problem
[nns,m] = size(A);
% Initialization of the arrays
J=zeros(nns,3);
V=zeros(nns,3);
% Definning for all the rows of H, the indexes j different than 0 (jj) and
% the values of A at these indexes (vv) in the arrays J and V.
for i=1:nns
    [~,jjf,vvf]=find(A(i,:));
    [vvn,vvi]=sort(vvf,'descend');
    ll=min(length(jjf),3);
    vv=vvn(1:ll);
    jj=jjf(vvi(1:ll));
    [jjs,jji]=sort(jj);
    vvs=vv(jji);
    J(i,1:ll)=jjs;
    V(i,1:ll)=vvs;
end
% New row indexes
i=(1:m)';
% Find the maximum of the non-zero values of each row
[yy,xx,~]=find((V==max(V,[],2))');
ind=sub2ind(size(J),xx,yy);
% New column indexes
j=J(ind);
% New values
v=ones(size(i));
% New sparse matrix
H=sparse(i,j,v,m,n);

end

