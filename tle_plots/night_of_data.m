function night_of_data(night,j)


cameras=['A' 'B' 'C' 'D' 'E' 'F' 'G'];
colors=['r', 'b', 'g', 'k','m', 'y', 'c'];

for i=1:length(cameras)
    AA=textread(['f' cameras(i) night 'M.fix']);
    figure(j)
    plot((AA(:,8)),AA(:,9),['.' colors(i)])
    hold on     
end
hold on 
title(['Data on ' night])
xlabel('azimuth (Degs)');
ylabel('elevation (Degs)');
grid on;