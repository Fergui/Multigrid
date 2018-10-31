function sn=findmgfit(n)
% find subset index range that give multigrid levels
% find number n1<=n such that n1=c*2^k+1 with small c and n1 close to n
for c=5:12
    levels=floor(log((n-1)/c)/log(2));
    ns(c)=c*2^levels+1;
end
n_mg=max(ns);
sn=floor((n-n_mg)/2)+[1:n_mg];
