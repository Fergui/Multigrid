function kmlcell = kmls2cell(regex)
% From a regex of kml files, it creates a cell of structs

kmls=dir(regex);
nn=length(kmls);
kmlcell=cell(nn,1);
for i=1:nn
    kmlcell{i}=kml2struct(kmls(i).name);
    kmlcell{i}.Name=kmls(i).name;
    date=kmls(i).name(17:29);
    kmlcell{i}.Date=strcat(date(1:4),'-',date(5:6),'-',date(7:8),'T',date(10:11),':',date(12:13),':00-00:00');
end

end

