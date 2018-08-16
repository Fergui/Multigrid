function i=binsearch(s,T)
% return smallest i such that s(i)=T
% in
%   s   array sorted nondecreasing
%   T   target value
% out
%   i   smallest index i such that s(i)=T

l = 1;
u = length(s)
for k=1:ceil(log(u)/log(2))
    if s(l)  > T, i=[]; break; end
    if s(u)  < T, i=[]; break; end
    m = floor((l + u)/2);
    % fprintf('s(%g)=%g s(%g)=%g s(%g)=%g\n',l,s(l),m,s(m),u,s(u))
    if s(m)<T, 
        l=m+1; 
    elseif s(m)>T, 
        u=m-1;
    else
        i=m;
        break
    end
end