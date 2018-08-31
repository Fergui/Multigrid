%% Multigrid Ideal from WRF-SFIRE execution test
% It is required a wrfout ideal case file set in tests/data/wrfout

% Setup the case
[u,p] = setup_file('data/wrfout',50);

% Special configurations for the multigrid method
% Using penalty
p.penalty=1; 
% Using dynamic rate of spread 
p.ros=0; 
% Recording in a gif the optimization plots
p.rec=0; 
% Showing the plots
p.plt=1;
% Strategy vector for the multigrid method
maxs=4;
p.multigrid=zeros(1,maxs);
for k=1:size(p.multigrid,2)+1
    p.multigrid(k)=k;
end
p.multigrid=flip(p.multigrid);

s.u=u; s.p=p;
clear u p
% Run the multigrid method
[um,Jop] = multigrid_process(s);