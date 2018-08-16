s=[ 1     3     3     3     4     4     7     8     9    10];

for smin=0:0.5:11
    for smax=0:0.5:11
        ii=binsearch_bounds(s,smin,smax);
        jj=find(s<=smax & s>=smin);
        if length(ii)~=length(jj) || any(ii~=jj),
            disp('test failed')
        end
    end
end