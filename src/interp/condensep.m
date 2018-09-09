function [H,rows,prows] = condensep(A,b,T)
% Call:
% [H,rows] = condensep(A,b,T)
% [H,rows,prows] = condensep(A,b,T)
%
% Description:
% Condense the rows of H in a strategic way to have the right side compatible
%
% Inputs:
%   A     sparse matrix to condense
%   b     right hand side of constraints Ax=b
%   T     time we want in the conflicting places (both perimeters are at the same place)
% Outputs:
%   H     sparse matrix condensed
%   rows  rows of the initial matrix A of the final matrix H
%   prows final rows of A which stays in H and were conflicting places
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Removing rows with empty intersection of nodes depending on
% Initialization
l=length(find(A(1,:)));
[jj,ii,~]=find(A');
% Definning for all the rows of H, the indexes j different than 0 (jj) 
% in the array J.
J=(reshape(jj,l,length(jj)/l))';
% Definning the matrix Jb with the values of the right hand
Jb=[J,b];
% Defining a cell of arrays where each cell is going to have problematic
% rows
fr=cell(l,1);
fn=cell(l,1);
for k=1:l
    % Defining the submatrices with one column of node index and right hand side
    Jk=Jb(:,[k end]);
    % Sorting the matrix with respect to the first column
    [JS,ind]=sortrows(Jk);
    % Distance matrix between consecutive rows
    D=abs(JS(2:end,:)-JS(1:end-1,:));
    % Find all the places with change of right hand side
    ip=find(D(:,2));
    % Find all the places with change of right hand side and with the same
    % node number
    rp=unique(ip(D(ip,1)==0));
    % Find which ones we haven't the right hand side that we want to conserve
    kk=find(JS(rp,2)~=T);
    % Defining the rows from the initial J where there is a conflict and we
    % don't have the right hand side that we want
    r=rp(kk);
    rp(kk)=[];
    rr=[r;rp+1];
    rn=[r+1;rp];
    fr{k}=ind(rr);
    fn{k}=ind(rn);
end
% Taking the rows that we want to remove from H and g
rro=[];
nro=[];
for k=1:l
    rro=unique([rro;fr{k}]);
    nro=unique([nro;fn{k}]);
end
rows=ii;
rows(rro)=[];
% Removing the rows from H and g
A(rro,:)=[];

%% Only rows linearly independent
H=A;
ixh=indep(H');
nrro=setdiff(1:length(ii),rro)';
nixh=ixh;
nixh(setdiff(nro,nrro))=[];
prows=rows(nixh);
rows=rows(ixh);
ind=setdiff(1:size(H,1),ixh);
H(ind,:)=[];
end

