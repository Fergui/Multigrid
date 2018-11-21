function varargout = repre(A)
% Call:
% H = repre(A)
%
% Description:
% Take one unique non-zero value (representative) for each row of A
%
% Call:
% [H,rows] = repre(A)
%
% Description:
% Take one unique non-zero value (representative) for each row of A and
% compute the ones which are linearly independent.
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
[m,n] = size(A);
% Initialization of the arrays
J=zeros(m,3);
V=zeros(m,3);
% Definning for all the rows of H, the indexes j different than 0 (jj) and
% the values of A at these indexes (vv) in the arrays J and V.
for i=1:m
    [~,jj,vv]=find(A(i,:));
    ll=length(jj);
    J(i,1:ll)=jj;
    V(i,1:ll)=vv;
end
% New row indexes
i=(1:m)';
% Find the maximum of the non-zero values of each row
[yy,xx,~]=find((V==max(V,[],2))');
kk=find(xx(2:end)-xx(1:end-1)==0);
xx(kk)=[];
yy(kk)=[];
ind=sub2ind(size(J),xx,yy);
% New column indexes
j=J(ind);
% New values
v=ones(size(i));
% New sparse matrix
H=sparse(i,j,v,m,n);

%% Only rows linearly independent
if nargout>1
    ixh=indep(H');
    rows=i(ixh);
    ind=setdiff(1:size(H,1),ixh);
    H(ind,:)=[];
    varargout{2}=rows;
end
varargout{1}=H;

end

