%% Test setup_real.m function

[u,p,s] = setup_real('./data/in.mat');

pr.u=u;
clear u
pr.p=p;
clear p
pr.s=s;
clear s

save('setup.mat','-struct','pr','-v7.3');

