%% Setup
mi=100; ni=100; dx=1; dy=1;
p=setup_circles_cascadic(mi,ni,dx,dy);
times=[p.per1_time,p.per2_time];

%% Assembly of the Cascadic Multigrid
[mg,levels]=assembly_cascadic_multigrid(p,times);

%% Strategy vector
exp=0:levels-1;
strategy=5.^exp;

%% Cascadic Multigrid
mg=cascadic_multigrid(mg,strategy);

%% Visualization of the result
U=mg(1).u; U(~mg(1).vmask)=nan;
figure, mesh(U), view([0,1]), title('Final result after the cascadic multigrid');
