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

% Computing lenght of non-zero values of each row of A
l=length(find(A(1,:)));
% Dimensions of A
[m,n]=size(A);
% Find non-zero values of A
[ii,jj,vv]=find(A);
% Sort the non-zero values by rows
[~,rows]=sort(ii);
% Sort the columns and the values
jjs=jj(rows);
vvs=vv(rows);
% Generate the matrix of columns and of values
J=(reshape(jjs,l,length(jjs)/l))';
V=(reshape(vvs,l,length(vvs)/l))';
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

