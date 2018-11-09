function C=coarsening(A,ratio)

[m,n]=size(A);
AC=zeros(m,n);
% Mean of the interior of A
Ac=A(2:end-1,2:end-1);
Aw=A(1:end-2,2:end-1);
Ae=A(3:end,2:end-1);
An=A(2:end-1,3:end);
As=A(2:end-1,1:end-2);
Asw=A(1:end-2,1:end-2);
Ase=A(3:end,1:end-2);
Ane=A(3:end,3:end);
Anw=A(1:end-2,3:end);
AM=(Ac+Aw+Ae+An+As+Asw+Ase+Ane+Anw)/9;
AC(2:end-1,2:end-1)=AM;
% Mean of the sud edge
Sc=A(2:end-1,1);
Sw=A(1:end-2,1);
Se=A(3:end,1);
Sn=A(2:end-1,2);
Sne=A(3:end,2);
Snw=A(1:end-2,2);
SM=(Sc+Sw+Se+Sn+Sne+Snw)/6;
AC(2:end-1,1)=SM;
% Mean of the nord edge
Nc=A(2:end-1,end);
Nw=A(1:end-2,end);
Ne=A(3:end,end);
Ns=A(2:end-1,end-1);
Nse=A(3:end,end-1);
Nsw=A(1:end-2,end-1);
NM=(Nc+Nw+Ne+Ns+Nse+Nsw)/6;
AC(2:end-1,end)=NM;
% Mean of the weast edge
Wc=A(1,2:end-1);
Ws=A(1,1:end-2);
Wn=A(1,3:end);
We=A(2,2:end-1);
Wse=A(2,1:end-2);
Wne=A(2,3:end);
WM=(Wc+Ws+Wn+We+Wse+Wne)/6;
AC(1,2:end-1)=WM;
% Mean of the east edge
Ec=A(end,2:end-1);
Es=A(end,1:end-2);
En=A(end,3:end);
Ew=A(end-1,2:end-1);
Esw=A(end-1,1:end-2);
Enw=A(end-1,3:end);
EM=(Ec+Es+En+Ew+Esw+Enw)/6;
AC(end,2:end-1)=EM;
% Mean of sud-weast corner
AC(1,1)=(A(1,1)+A(2,1)+A(2,2)+A(1,2))/4;
% Mean of sud-east corner
AC(end,1)=(A(end,1)+A(end-1,1)+A(end-1,2)+A(end,2))/4;
% Mean of nord-weast corner
AC(1,end)=(A(1,end)+A(1,end-1)+A(2,end-1)+A(2,end))/4;
% Mean of nord-east corner
AC(end,end)=(A(end,end)+A(end-1,end)+A(end-1,end-1)+A(end,end-1))/4;
% Final coarsened matrix
C=AC(1:ratio:m,1:ratio:n);

end