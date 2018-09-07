function kmlcell = kmls2cell(regex)
% Call:
% kmlcell = kmls2cell(regex)
%
% Description:
% From a regex of kml files, it generates an array of structs including all the kml files.
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
kmlcell=kml2struct(strcat(kmls(1).folder,'/',kmls(1).name));
for i=2:nn
    file=strcat(kmls(i).folder,'/',kmls(i).name);
    kml=kml2struct(file);
    kmlcell=[kmlcell;kml];
end
    
end

