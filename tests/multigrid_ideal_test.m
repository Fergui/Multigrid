%% Multigrid Ideal test

% Setup the case
[u,p] = setup_ideal(100,100,1,1,@s_eno,'c');
s.u=u; s.p=p;
clear u p

% Run the multigrid method
[um,Jop] = multigrid_process(s);