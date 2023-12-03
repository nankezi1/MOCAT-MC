function [datenumber] = date2num(date_str)
%DATE2NUM 
%   Converts date string to fractional year

% date_str = '2017-12-26T18:57:36+00:00'; % datestr(___,'yyyy-mm-ddTHH:MM:SS');

ndates = size(date_str,1);
datenumber = zeros(ndates,1);

for ind = 1:ndates
    curds = date_str(ind,:);
    year = str2double(curds(1:4));
    month = str2double(curds(6:7));
    day = str2double(curds(9:10));
    hour = str2double(curds(12:13));
    minute = str2double(curds(15:16));
    second = str2double(curds(18:19));
    JD = juliandate(year,month,day,hour,minute,second);
    JD_end = juliandate(year+1,01,01,00,00,00);
    JD_start = juliandate(year,01,01,00,00,00);

    datenumber(ind,1) = year+(JD-JD_start)/(JD_end-JD_start);
end

end

