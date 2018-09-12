%% Multigrid Real from WRF-SFIRE execution test

%% Case parameters
sim_path='/uufs/chpc.utah.edu/common/home/kochanski-group1/AIRPACT/Cougar_Creek';
per_path='/uufs/chpc.utah.edu/common/home/kochanski-group1/AIRPACT/Cougar_Creek_Perimeters/kml_files/wa*_dd83.kml';
itime='2015-08-11_01:00:00';
ilon=-121.374;
ilat=46.134;
dyn=0;

%% Preprocessing the data
if exist('data/setup.mat','file')
    fprintf('Loading data/setup.mat file...\n');
    load data/setup.mat
else
    if exist('data/in.mat','file')
        copyfile data/in.mat ./in.mat
    else
        fprintf('Preprocessing the data...\n');
        preproc(sim_path,per_path,itime,ilon,ilat,dyn);
    end 
    fprintf('Setting up the case...\n'); 
    [u,p,s]=setup_real('in.mat');
    pr.u=u; pr.p=p; pr.s=s;
    save('setup.mat','-struct','pr','-v7.3');
    clear pr
end

%% Experiment parameters
% Dynamic rate of spread?
p.ros=0; 
% Record in a gif file?
p.rec=0; 
% Show the plots?
p.plt=0;

%% Special configurations for the multigrid method
p.exp='real';
p.max_step=1.0;
p.nmesh=5;
p.max_depth=2;
p.min_depth=1;
p.mcycle=4;
p.penalty=1; 
maxs=8;
p.multigrid=zeros(1,maxs);
for k=1:size(p.multigrid,2)+1
    p.multigrid(k)=1;
end
p.multigrid=flip(p.multigrid);

%% Running the multigrid method
pr=s;
clear s
s.u=u; s.p=p; s.s=pr;
clear u p pr
% Run the multigrid method
[um,Jop] = multigrid_process(s);

%% Save results
pr.um=um; pr.Jop=Jop;
save('out.mat','-struct','pr','-v7.3');
fprintf('SUCCESS\n');
