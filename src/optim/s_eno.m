function [diff2,diff2du] = s_eno(diffL,diffR,ds)
% diff2 = s_eno(diffL,diffR)
% [diff2,diff2du] = s_eno(diffL,diffR,ds)
%input:
%   diffL       left difference ((u(i)-u(i-1))/ds)
%   diffR       right difference ((u(i+1)-u(i))/ds)
%   ds          distance between nodes
%ouput:
%   diff2    derivate that ENO method computes
%   diff2du  structure, with:
%           lp  one sided derivative from the positive side of diff2 respect to ui-1 (left node)
%           ln  one sided derivative from the negative side of diff2 respect to ui-1 (left node)
%           cp  one sided derivative from the positive side of diff2 respect to ui (center node)
%           cn  one sided derivative from the negative side of diff2 respect to ui (center node)
%           rp  one sided derivative from the positive side of diff2 respect to ui+1 (right node)
%           rn  one sided derivative from the negative side of diff2 respect to ui+1 (right node)

[m,n]=size(diffL);
diff2=zeros(m,n);
z=zeros(m,n);
diff2du=struct('lp',z,'ln',z,'cp',z,'cn',z,'rp',z,'rn',z);
C=(diffL<0.).*(diffR<0.);
diff2=diff2+C.*diffR;
diff2du.cp=diff2du.cp+C*(-1/ds);
diff2du.cn=diff2du.cn+C*(-1/ds);
diff2du.rp=diff2du.rp+C*(1/ds);
diff2du.rn=diff2du.rn+C*(1/ds);
C=(diffL>0.).*(diffR>0.);
diff2=diff2+C.*diffL;
diff2du.lp=diff2du.lp+C*(-1/ds);
diff2du.ln=diff2du.ln+C*(-1/ds);
diff2du.cp=diff2du.cp+C*(1/ds);
diff2du.cn=diff2du.cn+C*(1/ds);
C=(diffL>=0.).*(diffR<0.).*(diffL<-diffR); % abs(diffR)>=abs(diffL)
diff2=diff2+C.*diffR;
diff2du.cp=diff2du.cp+C*(-1/ds);
diff2du.cn=diff2du.cn+C*(-1/ds);
diff2du.rp=diff2du.rp+C*(1/ds);
diff2du.rn=diff2du.rn+C*(1/ds);    
C=(diffL>0.).*(diffR<=0.).*(diffL>-diffR);
diff2=diff2+C.*diffL;
diff2du.lp=diff2du.lp+C*(-1/ds);
diff2du.ln=diff2du.ln+C*(-1/ds);
diff2du.cp=diff2du.cp+C*(1/ds);
diff2du.cn=diff2du.cn+C*(1/ds);
C=(diffL>0.).*(diffR<0.).*(diffL==-diffR);
diff2=diff2+C.*diffL;
diff2du.ln=diff2du.ln+C*(-1/ds);
diff2du.cp=diff2du.cp+C*(1/ds);
diff2du.cn=diff2du.cn+C*(1/ds);
diff2du.rn=diff2du.rn+C*(1/ds);
C=(diffL<0.).*(diffR==0.);
diff2du.cp=diff2du.cp+C*(-1/ds);
diff2du.rn=diff2du.rn+C*(1/ds);
C=(diffL==0.).*(diffR>0.);
diff2du.ln=diff2du.ln+C*(-1/ds);
diff2du.cp=diff2du.cp+C*(1/ds);
C=(diffL==0.).*(diffR==0.);
diff2du.ln=diff2du.ln+C*(-1/ds);
diff2du.cp=diff2du.cp+C*(1/ds);
diff2du.rn=diff2du.rn+C*(1/ds);
end

