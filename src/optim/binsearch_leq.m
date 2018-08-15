function i=binsearch_leq(s,T)
% Call:
% i=binsearch_leq(s,T)
%
% Description:
% Return smallest i such that T <= s(i)
%
% Inputs:
%   s   array sorted nondecreasing
%   T   target value
% Outputs:
%   i   minimal such that T <= s(i)
%
%-------------------------------------------------------------------------

% j > min{ i: s(i) >=T} <=> exists i<j : s(i)>=T <=> (j > 1 && s(j-1) >= T)
% j < min{ i: s(i) >=T} <=> s(j) < T

l = 1;
i=0;   % initialize invalid value
u = length(s);
if T > s(u)
    i=[];
else
    for k=1:u
        % if s(l) > T, i=[]; break; end
        % not needed: if (l > 1 && s(l-1) >= T) i=[]; break; end
        % not needed: if s(u) < T, i=[]; break; end
        m = floor((l + u)/2);
        % fprintf('s(%g)=%g s(%g)=%g s(%g)=%g\n',l,s(l),m,s(m),u,s(u))
        if s(m) < T
            l=m+1; 
        % elseif s(m)>T, 
        elseif (m > 1 && s(m-1) >= T)
            u=m-1;
        else
            i=m;
            break
        end
    end
    if i==0
        error('binsearch_leq: bailed out of infinite loop - something wrong')
    end
    if ~ (T <= s(i) && (i==1 || T > s(i-1)))
        error('binsearch_leq: result failed output test')
    end
end
% test
end
