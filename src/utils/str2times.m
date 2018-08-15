function time=str2times(stri,strf)
% Call:
% time=str2time(stri,strf)

% Description:
% Compute the time in seconds from stri date to strf date

% Inputs:
%   stri    initial char string with YYYY-MM-DD_hh:mm:ss format
%   strf    final char string with YYYY-MM-DD_hh:mm:ss format

% Outputs:
%   time    difference of seconds between dates

% Developed in Matlab 9.2.0.556344 (R2017a) on MACINTOSH. 
% Angel Farguell (angel.farguell@gmail.com), 2018-08-15
%-------------------------------------------------------------------------

stri=char(stri);
strf=char(strf);
y1=str2double(stri(1:4));
y2=str2double(strf(1:4));
ms1=isleap(y1);
ms2=isleap(y2);
m1=str2num(stri(6:7));
m2=str2num(strf(6:7));
dd=str2double(strf(9:10))-str2double(stri(9:10));
dh=str2double(strf(12:13))-str2double(stri(12:13));
dmi=str2double(strf(15:16))-str2double(stri(15:16));
ds=str2double(strf(18:19))-str2double(stri(18:19));
dy=y2-y1;
dm=m2-m1;
if dy ~= 0 
    dd=sum(ms1(m1+1:end))+sum(ms2(1:m2-1))+ms1(m1)+dd;
end
if dm ~= 0
    dd=sum(ms1(m1+1:end))-sum(ms1(m2+1:end))+ms1(m1)-ms1(m2)+dd;
end
time=1440*dd+60*dh+dmi+ds/60;

    function m=isleap(y)
        if mod(y,4)==0 && mod(y,100)~=0 || mod(y,400)==0
            m=[31,29,31,30,31,30,31,31,30,31,30,31];
        else
            m=[31,28,31,30,31,30,31,31,30,31,30,31];
        end
    end
end
