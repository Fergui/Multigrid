function ii=binsearch_bounds(s,smin,smax)
% Call:
% ii=binsearch_bounds(s,smin,smax)
%
% Description:
% Search the array of all indices such that smin <= s(i) <= smax
%
% Inputs:
%   s           array sorted in increasing order
%   smin,smax   bounds to select values of xs
% Outputs:
%   ii          array of all indices such that smin <= s(i) <= smax
%
%-------------------------------------------------------------------------

imin=binsearch_leq(s,smin);
if isempty(imin),
    ii=[];
    return
end
imax=length(s)+1-binsearch_leq(flip(-s),-smax);
if isempty(imax),
    ii=[];
    return
end
ii=imin:imax;
end
