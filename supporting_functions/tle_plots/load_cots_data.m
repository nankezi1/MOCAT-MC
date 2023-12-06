close all; 
clear all; 
clc; 


night='140922';
night_of_data(night,1)

cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
cameras=['F']
plot_sat_data(28884,10)

% plot_sat_data(night,cameras,'GALAXY 15',11,10)


% % F 11 47 78, 56
% load('0922GL15')
% 
% for i=1:58
% cameras=['A'];
% plot_sat_data(night,cameras,'GALAXY 15',11,10)
% end

