%% Multigrid Ideal from WRF-SFIRE execution test
% It is required a wrfout ideal case file set in tests/data/wrfout

% Setup the case
[u,p] = setup_file('data/wrfout',50);
s.u=u; s.p=p;
clear u p

% Special configurations for the multigrid method
p.max_step=1.0;
p.nmesh=5;
p.min_depth=1;
p.max_depth=2;
% Number of multigrid cycles
p.mcycle=10; 
% Using penalty
p.penalty=1; 
% Using dynamic rate of spread 
p.ros=1; 
% Recording in a gif the optimization plots
p.rec=0; 
% Strategy vector for the multigrid method
maxs=4;
p.multigrid=zeros(1,maxs);
for k=1:size(p.multigrid,2)+1
    p.multigrid(k)=4;
end

% Run the multigrid method
[um,Jop] = multigrid_process(s);