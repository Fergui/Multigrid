%% Multigrid Ideal test of the exponent p.q
% Any data file is required

% Setup the case
[u,p] = setup_ideal(100,100,1,1,@s_eno,'c');

% Special configurations for the multigrid method
p.max_step=1.0;
p.nmesh=5;
p.min_depth=1;
p.max_depth=2;
% Number of multigrid cycles
p.mcycle=4; 
% Using penalty
p.penalty=1; 
% Using dynamic rate of spread 
p.ros=0; 
% Recording in a gif the optimization plots
p.rec=0; 
% Showing the plots
p.plt=0;
% Strategy vector for the multigrid method
maxs=3;
p.multigrid=zeros(1,maxs);
NN=size(p.multigrid,2)+1;
for k=1:NN
    p.multigrid(k)=NN-(k-1);
end
p.multigrid

% Experiment using different values of p.q
pr.u=u;
s=cell(1,5);
for kk=1:5
    p.q=2^kk;
    pr.p=p;
    s{kk}.p=p;
    % Run the multigrid method
    [s{kk}.um,s{kk}.Jop] = multigrid_process(pr); 
end
