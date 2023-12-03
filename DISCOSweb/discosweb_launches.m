% Gather LAUNCH data from DISCOS via discosweb: https://discosweb.esoc.esa.int/

% Created Sept 2022 - D.Jang
% May need updating (see discosweb.m)


% d = webread('https://discosweb.esoc.esa.int/api/objects')
%  https://discosweb.esoc.esa.int/apidocs/v2#section/Getting-Started

% Below is Dan's ESA API token - please use appropriately! 
token = 'ImI5MjBiMzM2LTk2ZjktNDA0ZC1hZGFiLTNjMmFiZjc4NGFkYSI.0_k20tZIVHaaq-89tngta3rmZd0';


%% curl version
% curl --header "Authorization: Bearer <token>"  --header "DiscosWeb-Api-Version: 2" "https://discosweb.esoc.esa.int/api/objects?filter=eq(objectClass,Payload)"
% Example: weboptions('HeaderFields',{'Content-Length' '78';'Content-Type' 'application/json'}) creates a weboptions object that contains two header fields: Content-Length with value 78 and Content-Type with value application/json.
wo = weboptions("HeaderFields",{'Authorization', ['Bearer ' token]; 'DiscosWeb-Api-Version' '2'});
initialdate = '2010-01-01'; %inclusive
finaldate = '2017-12-31'; %inclusive
filter_str = ['ge(epoch,epoch:''',initialdate,''')&le(epoch,epoch:''',finaldate,''')']; %greater than and less than specified dates
url_main = 'https://discosweb.esoc.esa.int/api/launches';
jsondata = webread(url_main,'filter',filter_str,'sort','epoch','page[size]','100','page[number]','1',wo);

% jsondata = webread(url_temp);
numpages = jsondata.meta.pagination.totalPages;

clear d;
launch_data = jsondata.data; %check d for structure fields

%% need to do 20 pages at a time, due to "TOO MANY REQUESTS" error 
for ind = 1:numpages
    if mod(ind,20)==1
        pause(60); %pause for 60 seconds before querying next 20 pages
    end

    jsondata = webread(url_main,'filter',filter_str,'sort','epoch','page[size]','100','page[number]',num2str(ind),wo);
    launch_data = [launch_data; jsondata.data];
    fprintf('%i ',ind);    
end
save('DISCOS_launches.mat','launch_data')

launchdates = zeros(length(launch_data),1);
for i=1:length(launch_data)
    launchdates(i) = date2num(launch_data(i).attributes.epoch);
end
figure
histogram(launchdates)

payload_discos = [118,127,128,204,193,175,175,392]; %from 2010 to end of 2017
average_payload_discos = sum(payload_discos)/numel(payload_discos);

launchespaper = [246.03174603174605, 0.7692307692307665;
298.41269841269843, 0.8371040723981906;
347.61904761904765, 2.058823529411768;
400, 5.1809954751131215;
450.7936507936508, 43.529411764705884;
500, 15.565610859728508;
550.7936507936508, 14.004524886877824;
600, 22.013574660633484;
650.7936507936508, 15.904977375565608;
700, 8.3710407239819;
749.2063492063492, 13.393665158371043;
800, 4.773755656108598;
850.7936507936507, 1.990950226244344;
900, 1.7194570135746616;
949.2063492063492, 1.312217194570131;
1001.5873015873016, 1.312217194570131;
1052.3809523809523, 1.312217194570131;
1100, 1.312217194570131;
1150.7936507936506, 1.7873303167420787;
1200, 2.058823529411768;
1250.7936507936506, 0.4977375565610913;
1300, 0.4977375565610913;
1349.2063492063494, 0.4977375565610913;
1400, 2.805429864253398;
1449.2063492063494, 4.027149321266968;
1500, 1.7873303167420787]; %https://apps.automeris.io/wpd/

num_launches_DAMAGE = sum(launchespaper(:,2));

% all saved on 8/27/2022

%% manipulate data (match via NORAD id etc)
% some don't have some fields, including satno, so clean up data files

clear all; 
load d_all.mat;         % discosweb data from above (Aug 2022)
load all_2022.mat;      % all_2022 tle (Aug 2022) - use for orbiting satsz

% manipulate discoweb mass/size data
a = [d.attributes];      % relevant fields: satno, objectClass, mass, shape, diameter, span, height, width
aSatids = nan(size(a));  % preallocate for sat ids in discoweb (maintain index)
for ind = 1:numel(a)
    if ~isempty(a(ind).satno)
        aSatids(ind) = a(ind).satno;
%     else
%         disp('empty')
    end
end

% manipulate TLE data (sats)
tles = [sats{:}];
validsatIDs = [tles.satnum];

% make db of: ID, mass, size(s)
for ind = 1:numel(tles)
    curID = tles(ind).satnum;

%     tic
%     fun = @(x) ~isempty(a(x).satno) && a(x).satno == curID;
%     tf2 = arrayfun(fun, 1:numel(a));
%     aind = find(tf2);
%     toc
    
%     tic
    if ~isempty(curID)
        aind = aSatids == curID;
    else
        aind = 0;
    end
%     toc

    if any(aind)
        tles(ind).objectClass = a(aind).objectClass;
        tles(ind).mass = a(aind).mass;
        tles(ind).shape = a(aind).shape;
        tles(ind).diameter = a(aind).diameter;
        tles(ind).span = a(aind).span;
        tles(ind).height = a(aind).height;
        tles(ind).width = a(aind).width;  
        tles(ind).depth = a(aind).depth; 
    else
        warning('Satellite ID %i not found in discoweb data (ind: %i); skipping... \n',curID,ind)
    end
    if mod(ind,round(numel(tles)/10)) == 0
        fprintf('%0.1f%% ',100*ind/numel(tles));
    end
end

% save('mergedTLEs.mat','tles')

%% Filter data out into relevant portions
% LEO orbits, payload types, histogram of sizes and masses etc

clear
tic
load mergedTLEs.mat
toc  % 6 sec

tic
maxLEOradius = 3000;  % limit for apogee of "LEO" [km]
fun = @(x) tles(x).alta < (maxLEOradius / earthRadius('km'));
LEOind = arrayfun(fun, 1:numel(tles));
tlesLEO = struct2table(tles(LEOind));       % table

% divide into object classes
token = 'Ijg1OGY0MzVmLWUxZDMtNGQ0NC1hZjdkLTc2YjkyM2U1YzNmZCI.gDp-NvNMQeTYzmp_no3vUnANiZ8';
wo = weboptions("HeaderFields",{'Authorization', ['Bearer ' token]; 'DiscosWeb-Api-Version' '2'});
json = webread('https://discosweb.esoc.esa.int/api/object-classes',wo);
tmp = [json.data.attributes];
objectClasses = {tmp.name};

tlesLEOclass = struct;
for ind = 1:numel(objectClasses)
    fun = @(x) strcmp(tlesLEO.objectClass{x}, objectClasses{ind});
    inds = arrayfun(fun, 1:size(tlesLEO,1));
    tlesLEOclass.(erase(objectClasses{ind}," ")) = tlesLEO(inds,:);
end
toc % 21 sec


% save('mergedTLEs_LEOclasses.mat','tlesLEOclass')


%% analyze size weight etc: "CURRENT" TLE FOR SITUATION NOW
% load mergedTLEs_LEOclasses.mat

% compile masses per class; note that empty values are discarded
mPL = cell2mat(tlesLEOclass.Payload.mass);
% mPD = cell2mat(tlesLEOclass.PayloadDebris.mass);  % doesn't exist
mPFD = cell2mat(tlesLEOclass.PayloadFragmentationDebris.mass);
mPMRO = cell2mat(tlesLEOclass.PayloadMissionRelatedObject.mass);
mRB = cell2mat(tlesLEOclass.RocketBody.mass);
mRD = cell2mat(tlesLEOclass.RocketDebris.mass);
mRMRO = cell2mat(tlesLEOclass.RocketMissionRelatedObject.mass);
% mUNK = cell2mat(tlesLEOclass.Unknown.mass);

% max([mRB; mPL; mPD; mPFD; mUNK; mRMRO; mRD; mRMRO])

figure(1); clf; 
histogram(mPL(mPL<5000));           % outlier: ISS: 450000 kg
title(sprintf('Payload (total: %i)',numel(mPL)));
pause; clf;
histogram(mPFD);            % outlier
title(sprintf('payload Fragmented debris (total: %i)',numel(mPFD)));
pause; clf;
histogram(mPMRO); 
title(sprintf('payload mission related obj (total: %i)',numel(mPMRO)));
pause; clf;
histogram(mRB);            % outlier?
title(sprintf('Rocket bodies (total: %i)',numel(mRB)));
pause; clf;
histogram(mRMRO);             % outlier?
title(sprintf('Rocket mission related obj (total: %i)',numel(mRMRO)));






