clear all, close all, clc
%Read all the possible types of objectClass

load('mergedTLEs.mat')

idx_all = 1;
all_class{1} = tles(1).objectClass;
for i=2:numel(tles)
    classi = tles(i).objectClass;
    sum_eq = numel(all_class);
    for j=1:numel(all_class)
        if isequal(classi,all_class{j})
            sum_eq = sum_eq-1;
        end
    end
    if sum_eq==numel(all_class)
        idx_all = idx_all+1;
        all_class{idx_all} = classi;
    end
end
all_class{:}

%     'Payload'
%     'Rocket Body'
%     'Payload Mission Related Object'
%     'Rocket Mission Related Object'
%     'Rocket Fragmentation Debris'
%     'Payload Fragmentation Debris'
%     'Payload Debris'
%     'Rocket Debris'
%     'Other Debris'
%     'Unknown'


%                         Unknown
%                    RocketDebris
%                   PayloadDebris
%      PayloadFragmentationDebris
%                         Payload
%     PayloadMissionRelatedObject
%       RocketFragmentationDebris
%                      RocketBody
%      RocketMissionRelatedObject
%       OtherMissionRelatedObject
%                     OtherDebris