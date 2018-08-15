function [H,rows] = condense(A)
%[H,rows] = condense(A)
% Condense the rows of A to be one constraint for triangle and all linearly independent
% in
%   A     sparse matrix to condence
% out
%   H     sparse matrix condenced
%   rows  A indexes of final rows of H condenced

%% Mean of the rows in the same triangle
% Dimensions of the problem
l=length(find(A(1,:)));
[~,mm]=size(A);
[ii,jj,vv]=find(A);
[~,rows]=sort(ii);
jjs=jj(rows);
vvs=vv(rows);
% Definning for all the rows of H, the indexes j different than 0 (jj) and
% the values of A at these indexes (vv) in the arrays J and V.
J=(reshape(jjs,l,length(jjs)/l))';
V=(reshape(vvs,l,length(vvs)/l))';
% Sorting the rows depending on the j indexes different than 0 in A.
[JS,index]=sortrows(J);
[nns,~]=size(JS);
% Sorting the values.
VS=V(index,:);
% Computing the distance between the consecutive rows of JS.
D=abs(JS(2:end,:)-JS(1:end-1,:));
% Define the indexes where there is a change, starting with the first.
inn=find([1;sum(D,2)]);
inn=[inn;nns+1];
% Length of this number of indexes
ninn=length(inn);
% Initializating the arrays of indexes to define a new sparse matrix
rows=zeros(1,ninn-1);
isn=zeros(1,l*(ninn-1));
jsn=zeros(1,l*(ninn-1));
vsn=zeros(1,l*(ninn-1));
for i=1:ninn-1
    rows(i)=index(inn(i));
    isn(l*i-(l-1):l*i)=i*ones(1,l);
    jsn(l*i-(l-1):l*i)=JS(inn(i),:);
    vsn(l*i-(l-1):l*i)=mean(VS(inn(i):inn(i+1)-1,:),1);
end
k=find(jsn);
C=sparse(isn(k),jsn(k),vsn(k),ninn-1,mm);

%% Only rows linearly independent
H=C;
ixh=indep(H');
rows=rows(ixh);
ind=setdiff(1:size(H,1),ixh);
H(ind,:)=[];
end

