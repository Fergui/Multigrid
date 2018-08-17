function preproc(path,kmlfile,ilon,ilat)
% Call:
% preproc(path,kmlfile,fuels)
%
% Description:
% Preprocess all the important data from wrfout and kml to a matlab file
%
% Example: 
% preproc('/uufs/chpc.utah.edu/common/home/kochanski-group2/NASA_NOAA/Las_Conchas_4d_large/wrfout_d04*','doc.kml',-106.5421,35.8121);
%
% Inputs
%	path       path to simulation outputs wrfout files.
% 	kmlfile    kml file with the shapefiles of the perimeters.
%   ilon,ilat  lon,lat ignition coordinates
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

%% Compute times and ignition frame and time
wrfouts=dir(path);
nf=size(wrfouts,1);
if ~nf
   error('Error, the path to the wrfout files is incorrect.') 
end
nignited=1;
for i=1:nf
	wname=strcat(wrfouts(i).folder,'/',wrfouts(i).name);
    	if nignited
        	s=nc2struct(wname,{'Times','FIRE_AREA'},{});
        	if i==1
            		time=char(s.times)';
            		NN=size(time,1);
	    		nn=size(time,1);
            		p.sdates=cell(1,nf*nn);
            		p.stimes=zeros(1,nf*nn);
            		for k=1:nn
                		p.sdates{k}=time(k,:);
            		end
            		stri=p.sdates{1};
        	else
            time=char(s.times)';
            nn=size(time,1);
        end
        for t=1:nn
            if max(max(s.fire_area(:,:,t))) > 0
                ignF=(i-1)*NN+t;
                nignited=0;
                break;
            end
        end
    else
        s=nc2struct(wname,{'Times'},{});
        time=char(s.times)';
        nn=size(time,1);
    end
    for k=NN*(i-1)+1:NN*(i-1)+nn
    	p.sdates{k}=time(k-NN*(i-1),:);
    end
    for j=1:nn
    	p.stimes(NN*(i-1)+j)=str2times(stri,time(j,:));
    end
    if nn<NN
        p.sdates(NN*(i-1)+nn+1:end)=[];
        p.stimes=p.stimes(1:NN*(i-1)+nn);
    end
end
strs=p.sdates{ignF};
p.tig=str2times(stri,strs);

%% Take shapefiles information from kmlfile
kml=kml2struct(kmlfile);
p.kml=kml;
nkml=size(kml,2);
% Compute times from shapefiles (perimeters)
p.pdates=cell(1,nkml);
p.ptimes=zeros(1,nkml);
p.sframes=zeros(1,nkml);
for i=1:nkml
	p.pdates{i}=kml(i).Date;
 	p.ptimes(i)=str2times(stri,p.pdates{i});
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
wname=strcat(wrfouts(igfi).folder,'/',wrfouts(igfi).name);
p.ignS=nc2struct(wname,{'Times','FXLONG','XLONG','FXLAT','XLAT','FIRE_AREA','NFUEL_CAT','FMC_G','UF','VF','DZDXF','DZDYF'},{'DX','DY'},igfr);
p.ignS.coordinates=[ilon,ilat];

%% Take the necessary data at perimeter frames
for i=1:nkml
	if ~isnan(p.sframes(i))
		pfi=ceil(p.sframes(i)/NN); % In which file is the perimeter frame?
		pfr=p.sframes(i)-floor(p.sframes(i)/NN)*NN; % In which frame?
		% Some later time after the perimeter
		p.lf=20; % number of frames later
		framel=p.sframes(i)+p.lf;
		if framel>(nf-1)*NN+nn
			framel=(nf-1)*NN+nn;
		end
		pfil=ceil(framel/NN);
		pfrl=framel-floor(framel/NN)*NN;
		wname=strcat(wrfouts(pfi).folder,'/',wrfouts(pfi).name);
		wnamel=strcat(path,wrfouts(pfil).name);
		p.perS(i)=nc2struct(wname,{'Times','ROS','FIRE_AREA','FMC_G','UF','VF'},{},pfr);
		p.perlS(i)=nc2struct(wnamel,{'Times','ROS','TIGN_G','FIRE_AREA','FMC_G','UF','VF'},{},pfrl);
	end
end

%% Takin the dynamic variables at all the time steps
for i=1:nf
    wname=strcat(wrfouts(i).folder,'/',wrfouts(i).name);
    s=nc2struct(wname,{'UF','VF','FMC_G'},{});
    

%% Saving the final structure
save('setup.mat','p');
end

