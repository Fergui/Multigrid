function [H,g] = condense_arti(A,b,u)
% Call:
% [H,g] = condense_arti(A,b,u)
%
% Description:
% Condense the rows of H and creates an artificial perimeter
%
% Inputs:
%   A     sparse matrix to condense
%   b     right hand side of constraints Ax=b
%   u     fire arrival time to know where the solution is going
% Outputs:
%   H     sparse matrix condensed
%   g     right hand side compatible with the new H, s.t. Hx=g iff Ax=b
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

% Initialization
uu=unique(b);
l=length(find(A(1,:)));
[mm,nn]=size(A);
[jj,~,~]=find(A');
% Definning for all the rows of H, the indexes j different than 0 (jj) 
% in the array J.
J=(reshape(jj,l,length(jj)/l))';
% Definning the matrix Jb with the values of the right hand
Jb=[J,b];
% Sorting the matrix with respect to the first column
[JS,index]=sortrows(Jb);
% Distance matrix between consecutive rows
D=abs(JS(2:end,1:l)-JS(1:end-1,1:l));
S=sum(D,2);
% Find all the places with the same node representants
ip=find(S==0);
ip=[ip;ip+1];
% The rest of the rows (they are going to be the first rows of the new
% matrix H)
ir=(1:mm)';
ir(ip)=[];
rows1=index(ir);
% Generate the first part of the new constraint matrix
A1=A(rows1,:);
A1R=representant(A1);
[A1,r1]=condense(A1R);
rs1=rows1(r1);
b1=b(rs1);
% Compute the new part that we are going to re-define
Jn=JS(ip,:);
[JnS,kk]=sortrows(Jn);
rr2=index(ip(kk));
U=u(JnS(:,1:3));
[Umi,miJ]=min(U,[],2);
[Uma,maJ]=max(U,[],2);
miind = sub2ind(size(JnS), (1:size(JnS,1))', miJ);
maind = sub2ind(size(JnS), (1:size(JnS,1))', maJ);
Jmi=JnS(miind);
Jma=JnS(maind);
JnS=[Jmi,Jma,JnS(:,4)];
U=[Umi,Uma];
DU=U(:,2)-U(:,1);
% Look when u1<u2
ind=find(DU>0);
ind1=(JnS(ind,3)==uu(1));
ind2=(JnS(ind,3)==uu(2));
jj1=JnS(ind(ind1),1);
jj2=JnS(ind(ind2),2);
rs21=rr2(ind(ind1));
rs22=rr2(ind(ind2));
% Look when u1>u2
ind=find(DU<0);
ind1=(JnS(ind,3)==uu(1));
ind2=(JnS(ind,3)==uu(2));
jj1=[jj1;JnS(ind(ind1),2)];
jj2=[jj2;JnS(ind(ind2),1)];
rs21=[rs21;rr2(ind(ind1))];
rs22=[rs22;rr2(ind(ind2))];
% Look when u1==u2
ind=find(DU==0);
ind1=(JnS(ind,3)==uu(1));
jj1=[jj1;JnS(ind(ind1),1)];
rs21=[rs21;rr2(ind(ind1))];
% Generate A2 the second part of the new matrix
m21=length(jj1);
m22=length(jj2);
bb21=uu(1)*ones(m21,1);
bb22=uu(2)*ones(m22,1);
bb2=[bb21;bb22];
m2=m21+m22;
iis=(1:m2)';
jjs=[jj1;jj2];
rows2=[rs21;rs22];
vvs=ones(size(iis));
A2=sparse(iis,jjs,vvs,m2,nn);
if ~isempty(A2)
    [A2,r2]=condense(A2);
    rs2=rows2(r2);
    b2=bb2(r2);
    % Define final matrices
    An=[A1;A2];
    bn=[b1;b2];
    rows=[rs1;rs2];
else
    An=A1;
    bn=b1;
    rows=rs1;
end
% Condense final solution
[A,kk]=condensep(An,bn,uu(1));
b=bn(kk);
rows=rows(kk);


%% Only rows linearly independent
H=A;
g=b;
ixh=indep(H');
rows=rows(ixh);
ind=setdiff(1:size(H,1),ixh);
H(ind,:)=[];
g(ind)=[];
end

