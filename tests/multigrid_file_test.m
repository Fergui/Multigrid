%% Multigrid Ideal from WRF-SFIRE execution test
% It is required a wrfout ideal case file set in tests/data/wrfout

% Setup the case
[u,p] = setup_file('data/wrfout',50);
s.u=u; s.p=p;
clear u p

% Run the multigrid method
[um,Jop] = multigrid_process(s);