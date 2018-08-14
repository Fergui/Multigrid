function nfuelcat=fuels2fuel(nfuelscat)
% nfuelcat=fuels2fuel(nfuelscat)
% From matrix nfuelscat to struct of matrices with fuels type equal in A
% inputs:
%   nfuelscat      matrix with all the fuel types on it
% outputs:
%   nfuelcat       structure with a matrix in each fuel variable important
%                  to compute the ROS

run fuels
[m,n]=size(nfuelscat);
z=zeros(m,n);
nfuelcat=struct('windrf',z,'fgi',z,'fueldepthm',z,'savr',z,'fuelmce',z,'fueldens',z,'st',z,'se',z,'weight',z,'fci_d',z,'fct',z,'ichap',z,'fci',z,'fcbr',z,'hfgl',z,'cmbcnst',z,'fuelheat',z,'fuelmc_c',z);
for i=1:m
    for j=1:n
        if nfuelscat(i,j)==14
            % All the variables equal to 1 to not get Inf values
            nfuelcat.windrf(i,j)=1;
            nfuelcat.fgi(i,j)=1;
            nfuelcat.fueldepthm(i,j)=1;
            nfuelcat.savr(i,j)=1;
            nfuelcat.fuelmce(i,j)=1;
            nfuelcat.fueldens(i,j)=1;
            nfuelcat.st(i,j)=1;
            nfuelcat.se(i,j)=1;
            nfuelcat.weight(i,j)=1;
            nfuelcat.fci_d(i,j)=1;
            nfuelcat.fct(i,j)=1;
            nfuelcat.ichap(i,j)=0; % not chaparral because no fire
            nfuelcat.fci(i,j)=1;
            nfuelcat.fcbr(i,j)=1;
            nfuelcat.hfgl(i,j)=1;
            nfuelcat.cmbcnst(i,j)=0; % only variable needed to be 0 for ROS=0
            nfuelcat.fuelheat(i,j)=1;
            nfuelcat.fuelmc_c(i,j)=1;
        else
            s=fuel(nfuelscat(i,j));
            nfuelcat.windrf(i,j)=s.windrf;
            nfuelcat.fgi(i,j)=s.fgi;
            nfuelcat.fueldepthm(i,j)=s.fueldepthm;
            nfuelcat.savr(i,j)=s.savr;
            nfuelcat.fuelmce(i,j)=s.fuelmce;
            nfuelcat.fueldens(i,j)=s.fueldens;
            nfuelcat.st(i,j)=s.st;
            nfuelcat.se(i,j)=s.se;
            nfuelcat.weight(i,j)=s.weight;
            nfuelcat.fci_d(i,j)=s.fci_d;
            nfuelcat.fct(i,j)=s.fct;
            nfuelcat.ichap(i,j)=s.ichap;
            nfuelcat.fci(i,j)=s.fci;
            nfuelcat.fcbr(i,j)=s.fcbr;
            nfuelcat.hfgl(i,j)=s.hfgl;
            nfuelcat.cmbcnst(i,j)=s.cmbcnst;
            nfuelcat.fuelheat(i,j)=s.fuelheat;
            nfuelcat.fuelmc_c(i,j)=s.fuelmc_c;
        end
    end
end
end
        
        