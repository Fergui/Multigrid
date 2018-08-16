function [s,c]=getstartcount(varinfo)
% [s,c]=getstartcount(varinfo)
%
% Jan Mandel
%-------------------------------------------------------------------------

s=zeros(varinfo.ndims,1);
c=s;
for i=1:varinfo.ndims
   c(i)=varinfo.dimlength(i); 
end
end
