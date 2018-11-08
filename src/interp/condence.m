function [H,rows] = condence(A)
% Call:
% [H,rows] = condence(A)
%
% Description:
% Condence the rows of H to be one constraint for triangle and all linearly independent. 
% Slow but it works in all the cases.
%
% Inputs:
%   A     sparse matrix to condence
% Outputs:
%   H     sparse matrix condenced
%   rows  A indexes of final rows of H condenced
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Mean of the rows in the same triangle
% Dimensions of the problem
[nns,mm] = size(A);
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
% Sorting the rows depending on the j indexes different than 0 in A.
[JS,index]=sortrows(J);
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
isn=zeros(1,3*(ninn-1));
jsn=zeros(1,3*(ninn-1));
vsn=zeros(1,3*(ninn-1));
for i=1:ninn-1
    rows(i)=index(inn(i));
    isn(3*i-2:3*i)=i*ones(1,3);
    jsn(3*i-2:3*i)=JS(inn(i),:);
    vsn(3*i-2:3*i)=mean(VS(inn(i):inn(i+1)-1,:),1);
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

