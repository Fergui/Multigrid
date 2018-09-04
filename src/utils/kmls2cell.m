function kmlcell = kmls2cell(regex)
% Call:
% kmlcell = kmls2cell(regex)
%
% Description:
% From a regex of kml files, it creates a cell of structs including all the kml files.
%
% Inputs:
%   regex     path to a kml files using a regex expresion
% Outputs:
%   kmlcell   kml cell of structs for each kml
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-09-04
%-------------------------------------------------------------------------

kmls=dir(regex);
nn=length(kmls);
kmlcell=cell(nn,1);
for i=1:nn
    kmlcell{i}=kml2struct(kmls(i).folder+'/'+kmls(i).name);
    kmlcell{i}.Name=kmls(i).name;
    date=kmls(i).name(17:29);
    kmlcell{i}.Date=strcat(date(1:4),'-',date(5:6),'-',date(7:8),'T',date(10:11),':',date(12:13),':00-00:00');
end

end

