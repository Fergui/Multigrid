function kmlStruct = kml2struct(kmlFile)
% Call:
% kmlStruct = kml2struct(kmlFile)
%
% Description:
% Import a .kml file as a vector array of shapefile structs, with Geometry, Name, Description, Lon, Lat, and BoundaryBox fields.  
% Structs may contain a mix of points, lines, and polygons.
% .kml files with folder structure will not be presented as such, but will appear as a single vector array of structs.
%
% Inputs:
%       kmlFile    kml file
% Outputs:
%       kmlStruct  Structure with all the information inside the kml
%
% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
% Modified from above
%-------------------------------------------------------------------------

%Downloaded from https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct
%October 13, 2016
%
%Copyright (c) 2012, James Slegers
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%* Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in
%the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%                       SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%                       INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

[FID msg] = fopen(kmlFile,'rt');

if FID<0
    error(msg)
end

txt = fread(FID,'uint8=>char')';
fclose(FID);

expr = '<Placemark.+?>.+?</Placemark>';

objectStrings = regexp(txt,expr,'match');

Nos = length(objectStrings);

for ii = 1:Nos
    % Find Object Name Field
    bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
    if isempty(bucket)
        bucket = regexp(objectStrings{ii},'<SimpleData name="FIRE_NAME">.+?</SimpleData>','match');
        if isempty(bucket)
            name = 'undefined';
        else
            name = regexprep(bucket{1},'<SimpleData name="FIRE_NAME">\s*','');
            name = regexprep(name,'\s*</SimpleData>','');
        end
    else
        % Clip off flags
        name = regexprep(bucket{1},'<name.*?>\s*','');
        name = regexprep(name,'\s*</name>','');
    end
    
    % Find Object Date Field
    bucket = regexp(objectStrings{ii},'<TimeSpan.*?>.+?</TimeSpan>','match');
    if isempty(bucket)
        bucket = regexp(objectStrings{ii},'<SimpleData name="DATE_">.+?</SimpleData>','match');
        if isempty(bucket)
            date = 'undefined';
        else
            dat = regexprep(bucket{1},'<SimpleData name="DATE_">\s*','');
            dat = regexprep(dat,'\s*</SimpleData>','');
            bucket2 = regexp(objectStrings{ii},'<SimpleData name="TIME_">.+?</SimpleData>','match');
            time = regexprep(bucket2{1},'<SimpleData name="TIME_">\s*','');
            time = regexprep(time,'\s*</SimpleData>','');
            date = [dat(1:4),'-',dat(6:7),'-',dat(9:10),'T',time(1:2),':',time(3:4),':00-00:00'];
        end
    else
        % Clip off flags
        date = regexprep(bucket{1},'<TimeSpan.*?><begin>\s*','');
        date = regexprep(date,'\s*</begin></TimeSpan>','');
    end
    
    % Find Object Description Field
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    if isempty(bucket)
        desc = '';
    else
        % Clip off flags
        desc = regexprep(bucket{1},'<description.*?>\s*','');
        desc = regexprep(desc,'\s*</description>','');
    end
    
    geom = 0;
    % Identify Object Type
    if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
        geom = 1;
    elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
        geom = 2;
    elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
        geom = 3;
    end
    
    switch geom
        case 1
            geometry = 'Point';
        case 2
            geometry = 'Line';
        case 3
            geometry = 'Polygon';
        otherwise
            geometry = '';
    end
    
    % Find Coordinate Field
    bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
    for jj = 1:length(bucket)
       coordStr = regexprep(bucket{jj},'<coordinates.*?>(\s+)*','');
       coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
       coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
       [m,n] = size(coordMat);
       coordMat = reshape(coordMat,3,m*n/3)';
       bucket{jj} = coordMat;
    end

    % Create structure
    kmlStruct(ii).Geometry = geometry;
    kmlStruct(ii).Date = date;
    kmlStruct(ii).Name = name;
    kmlStruct(ii).Description = desc;
    kmlStruct(ii).Lon = cell(1,length(bucket));
    kmlStruct(ii).Lat = cell(1,length(bucket));
    kmlStruct(ii).BoundingBox = cell(1,length(bucket));
    for jj = 1:length(bucket)
        coord=bucket{jj};
        [Lat, Lon] = poly2ccw(coord(:,2),coord(:,1));
        if geom==3
            Lon = [Lon;NaN];
            Lat = [Lat;NaN];
        end
        kmlStruct(ii).Lon{jj} = Lon;
        kmlStruct(ii).Lat{jj} = Lat;
        kmlStruct(ii).BoundingBox{jj} = [[min(Lon) min(Lat);max(Lon) max(Lat)]];
    end
end
