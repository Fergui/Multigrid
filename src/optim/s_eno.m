function [varargout] = s_eno(diffL,diffR,ds)
% Call:
% diff2 = s_eno(diffL,diffR,ds)
%
% Description:
% Computes the ENO upwinding method.
%
% Call:
% [diff2,diff2du] = s_eno(diffL,diffR,ds)
%
% Description:
% Computes the ENO upwinding method and the second derivatives using the method proposed in Angel Farguell thesis.
%
% Inputs:
%   diffL       left difference ((u(i)-u(i-1))/ds)
%   diffR       right difference ((u(i+1)-u(i))/ds)
%   ds          distance between nodes
% Ouputs:
%   diff2    derivate that ENO method computes
%   diff2du  structure, with:
%           lp  one sided derivative from the positive side of diff2 respect to ui-1 (left node)
%           ln  one sided derivative from the negative side of diff2 respect to ui-1 (left node)
%           cp  one sided derivative from the positive side of diff2 respect to ui (center node)
%           cn  one sided derivative from the negative side of diff2 respect to ui (center node)
%           rp  one sided derivative from the positive side of diff2 respect to ui+1 (right node)
%           rn  one sided derivative from the negative side of diff2 respect to ui+1 (right node)
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

[m,n]=size(diffL);
diff2=zeros(m,n);
z=zeros(m,n);
C=(diffL<0.).*(diffR<0.);
diff2=diff2+C.*diffR;
C=(diffL>0.).*(diffR>0.);
diff2=diff2+C.*diffL;
C=(diffL>=0.).*(diffR<0.).*(diffL<-diffR); % abs(diffR)>=abs(diffL)
diff2=diff2+C.*diffR;  
C=(diffL>0.).*(diffR<=0.).*(diffL>-diffR);
diff2=diff2+C.*diffL;
C=(diffL>0.).*(diffR<0.).*(diffL==-diffR);
diff2=diff2+C.*diffL;
varargout{1}=diff2;

if nargout > 1
    diff2du=struct('lp',z,'ln',z,'cp',z,'cn',z,'rp',z,'rn',z);
    C=(diffL<0.).*(diffR<0.);
    diff2du.cp=diff2du.cp+C*(-1/ds);
    diff2du.cn=diff2du.cn+C*(-1/ds);
    diff2du.rp=diff2du.rp+C*(1/ds);
    diff2du.rn=diff2du.rn+C*(1/ds);
    C=(diffL>0.).*(diffR>0.);
    diff2du.lp=diff2du.lp+C*(-1/ds);
    diff2du.ln=diff2du.ln+C*(-1/ds);
    diff2du.cp=diff2du.cp+C*(1/ds);
    diff2du.cn=diff2du.cn+C*(1/ds);
    C=(diffL>=0.).*(diffR<0.).*(diffL<-diffR); % abs(diffR)>=abs(diffL)
    diff2du.cp=diff2du.cp+C*(-1/ds);
    diff2du.cn=diff2du.cn+C*(-1/ds);
    diff2du.rp=diff2du.rp+C*(1/ds);
    diff2du.rn=diff2du.rn+C*(1/ds);  
    C=(diffL>0.).*(diffR<=0.).*(diffL>-diffR);
    diff2du.lp=diff2du.lp+C*(-1/ds);
    diff2du.ln=diff2du.ln+C*(-1/ds);
    diff2du.cp=diff2du.cp+C*(1/ds);
    diff2du.cn=diff2du.cn+C*(1/ds);
    C=(diffL>0.).*(diffR<0.).*(diffL==-diffR);
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
    varargout{2}=diff2du;
end

end

