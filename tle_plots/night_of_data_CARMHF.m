function night_of_data_CARMHF(night,j)
cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
colors=['r', 'b', 'g', 'k','m', 'y', 'c'];
ym=[0.0662182180411231 -0.00132638607501781];
obs=[];

for i=1:length(cameras)
    AA=textread(['f' cameras(i) night 'M.fix']);
    obs=[obs;[90001*ones(length(AA(:,3)),1) AA(:,3)+2400000.5-2430000.0 (AA(:,6)-ym(1))*pi/180 (AA(:,7)-ym(2))*pi/180]];
    figure(j)
    plot((AA(:,8)),AA(:,9),['.' colors(i)])
    hold on     
end
% obs=obs(1:1000,:);
% dlmwrite('/Users/richlinares/Desktop/Rich/grants/2014/usno-cots/CAR-MHF/CARMHF_usno_v1/car-mhf-matlab/SampleObs/pseudo_GEODSS.gen',obs,'delimiter','\t','precision',10)
% 

hold on 
title(['Data on ' night])
xlabel('azimuth (Degs)');
ylabel('elevation (Degs)');
grid on;