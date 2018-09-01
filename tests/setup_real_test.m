%% Test setup_real.m function

% Run the setup_real function
fprintf('Running setup_real.m...\n');
[u,p,s] = setup_real('./data/in.mat');

% Save the results in a structure pr
fprintf('Saving the results...\n');
pr.u=u;
clear u
pr.p=p;
clear p
pr.s=s;
clear s
% Save in a file setup.mat
save('setup.mat','-struct','pr','-v7.3');
fprintf('SUCCESS\n');

