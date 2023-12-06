close all; 
clear all; 
clc; 


night='140922';
night_of_data_CARMHF(night,1)

 cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
% cameras=['F']
% plot_sat_data(28884,10)
catnum=[];
for i=1:length(cameras)
    [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10]=...
        textread(['xid.f' cameras(i) 'M'],...
        '%s %s %s %f %f %f %f %s %s %s', -1);
    catnum=[catnum;A5];
end
[catnum,ai,~]=unique(catnum);
% plot_sat_data(night,cameras,'GALAXY 15',11,10)
load_GEO_sats_CARMHF(catnum,night)

% % F 11 47 78, 56
% load('0922GL15')
% 
% for i=1:58
% cameras=['A'];
% plot_sat_data(night,cameras,'GALAXY 15',11,10)
% end

