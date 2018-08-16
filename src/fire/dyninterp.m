function ds=dyninterp(u,p)
% Temporal interpolation of the wind fields and the fuel moisture
% ds=dyninterp(u,p)
%in
%   u   fire arrival time
%   p   structure with:
%      vmask    matrix, true where the values of level set function can change
%      dynS     structure with: 
%           simt     array of simulation times in seconds
%           uf, vf   (x,y,t): wind components of the simulation
%           fmc_g    (x,y,t): fuel moisture of the simulation
%out
%   ds  dynam struct to compute the rate of spread, contains:
%       uf, vf   (x,y,t): wind coordinates interpolated
%       fmc_g    (x,y,t): fuel moisture interpolated

tt=p.dynS.simt;
[m,n]=size(u);
un=zeros(m,n);
vn=zeros(m,n);
fmcn=zeros(m,n);
for i=1:m
    for j=1:n
        if p.vmask(i,j)
            uij=squeeze(p.dynS.uf(i,j,:));
            un(i,j)=interp1(tt,uij,u(i,j));
            vij=squeeze(p.dynS.vf(i,j,:));
            vn(i,j)=interp1(tt,vij,u(i,j));
            fmcij=squeeze(p.dynS.fmc_g(i,j,:));
            fmcn(i,j)=interp1(tt,fmcij,u(i,j));
        end
    end
end

ds.uf=un;
ds.vf=vn;
ds.fmc_g=fmcn;

end