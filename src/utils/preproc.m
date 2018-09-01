function preproc(path,kmlfile,itime,ilon,ilat,dyn)
% Call:
% preproc(path,kmlfile,fuels)
%
% Description:
% Preprocess all the important data from wrfout and kml to a matlab file
%
% Example: 
% preproc('/glade/u/home/angelfc/project/lasconchas/simulation_large','doc.kml','2011-06-26_19:00:00',-106.542141,35.812056,1);
%
% Inputs
%	path       path to simulation outputs wrfout files.
% 	kmlfile    kml file with the shapefiles of the perimeters.
%	itime      ignition time
%   ilon,ilat  lon,lat ignition coordinates.
%   dyn        0: do not save dynS structure.
%              1: save dynS structure.
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Compute times and ignition frame and time
path=strcat(path,'/');
wrfouts=dir(strcat(path,'wrfout_d04*'));
nf=size(wrfouts,1);
nignited=1;
for i=1:nf
    wname=strcat(path,wrfouts(i).name);
    if dyn
        s=nc2struct(wname,{'Times','UF','VF','FMC_G'},{});
    else
        s=nc2struct(wname,{'Times'},{});
    end
    time=char(s.times)';
    nn=size(time,1);
    if i==1
        NN=size(time,1);
        kk=1;
        p.sdates=time;
        p.stimes=zeros(1,nf*NN);
        stri=time(1,:);
        if dyn
                p.dynS.u=s.uf(:,:,1:kk:end);
                s=rmfield(s,'uf');
                p.dynS.v=s.vf(:,:,1:kk:end);
                s=rmfield(s,'vf');
                p.dynS.fmc_g=s.fmc_g(:,:,1:kk:end);
                s=rmfield(s,'fmc_g');
        end
    end
    if nignited
        for t=1:nn
            if time(t,:)==itime
                ignF=(i-1)*NN+t;
                nignited=0;
                break;
            end
        end
    end
    if i>1
        p.sdates=[p.sdates;time];
        if dyn
                p.dynS.u=cat(3,p.dynS.u,s.uf(:,:,1:kk:end));
                s=rmfield(s,'uf');
                p.dynS.v=cat(3,p.dynS.v,s.vf(:,:,1:kk:end));
                s=rmfield(s,'vf');
                p.dynS.fmc_g=cat(3,p.dynS.fmc_g,s.fmc_g(:,:,1:kk:end));
                s=rmfield(s,'fmc_g');
        end
    end
    for j=1:nn
        p.stimes(NN*(i-1)+j)=str2time(stri,time(j,:));
    end
end
strs=p.sdates(ignF,:);
p.tig=str2time(stri,strs);

%% Take shapefiles information from kmlfile
kml=kml2struct(kmlfile);
dsize=size(kml(1).Date,2);
p.kml=kml;
nkml=size(kml,2);
% Compute times from shapefiles (perimeters)
p.pdates=cell(1,nkml);
p.ptimes=zeros(1,nkml);
p.sframes=zeros(1,nkml);
for i=1:nkml
        p.pdates{i}=kml(i).Date;
        p.ptimes(i)=str2time(stri,p.pdates{i});
        % Compute frames for perimeters
        post=find(p.stimes-p.ptimes(i)>0);
        if ~isempty(post)
                p.sframes(i)=post(1);
        else
                p.sframes(i)=nan;
        end
end

%% Take the necessary data at ignition point
igfi=ceil(ignF/NN); % In which file is the ignition frame?
igfr=ignF-floor(ignF/NN)*NN; % In which frame inside the file?
wname=strcat(path,wrfouts(igfi).name);
p.ignS=nc2struct(wname,{'Times','FXLONG','XLONG','FXLAT','XLAT','FIRE_AREA','NFUEL_CAT','FMC_G','UF','VF','DZDXF','DZDYF'},{'DX','DY'},igfr);
p.ignS.coordinates=[ilon,ilat];

%% Take the necessary data at perimeter frames
for i=1:nkml
        if ~isnan(p.sframes(i))
                pfi=ceil(p.sframes(i)/NN); % In which file is the perimeter frame?
                pfr=p.sframes(i)-floor(p.sframes(i)/NN)*NN; % In which frame?
                % Some later time after the perimeter
                lf=50; % number of frames later
                framel=p.sframes(i)+lf;
                if framel>(nf-1)*NN+nn
                        framel=(nf-1)*NN+nn;
                end
                pfil=ceil(framel/NN);
                pfrl=framel-floor(framel/NN)*NN;
                wname=strcat(path,wrfouts(pfi).name);
                wnamel=strcat(path,wrfouts(pfil).name);
                p.perS(i)=nc2struct(wname,{'Times','ROS','FIRE_AREA','FMC_G','UF','VF'},{},pfr);
                p.perlS(i)=nc2struct(wnamel,{'Times','ROS','TIGN_G','FIRE_AREA','FMC_G','UF','VF'},{},pfrl);
        end
end
    
%% Saving the final structure
save('setup.mat','p','-v7.3');
fprintf('SUCCESS\n');
end
