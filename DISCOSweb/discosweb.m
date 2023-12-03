% Gather OBJECT data from DISCOS via discosweb: https://discosweb.esoc.esa.int/
% Merge with TLE data from space-track.org
% Create initialized.mat in matsat format

% Created Sept 2022 - D.Jang
% Updated March 2023


%% 
% d = webread('https://discosweb.esoc.esa.int/api/objects')
%  https://discosweb.esoc.esa.int/apidocs/v2#section/Getting-Started

% Below is Dan's ESA API token - please use appropriately! 
token = 'ImI5MjBiMzM2LTk2ZjktNDA0ZC1hZGFiLTNjMmFiZjc4NGFkYSI.0_k20tZIVHaaq-89tngta3rmZd0';

%% python version
% from pprint import pprint
% import requests
% 
% URL = 'https://discosweb.esoc.esa.int'
% token = ''
% 
% response = requests.get(
%     f'{URL}/api/objects',
%     headers={
%         'Authorization': f'Bearer {token}',
%         'DiscosWeb-Api-Version': '2',
%     },
%     params={
%         'filter': "eq(objectClass,Payload)&gt(reentry.epoch,epoch:'2020-01-01')",
%         'sort': '-reentry.epoch',
%     },
% )
% 
% doc = response.json()
% if response.ok:
%     pprint(doc['data'])
% else:
%     pprint(doc['errors'])

%% curl version
% curl --header "Authorization: Bearer <token>"  --header "DiscosWeb-Api-Version: 2" "https://discosweb.esoc.esa.int/api/objects?filter=eq(objectClass,Payload)"
% Example: weboptions('HeaderFields',{'Content-Length' '78';'Content-Type' 'application/json'}) creates a weboptions object that contains two header fields: Content-Length with value 78 and Content-Type with value application/json.

wo = weboptions("HeaderFields",{'Authorization', ['Bearer ' token]; 'DiscosWeb-Api-Version' '2'});

jsondata = webread('https://discosweb.esoc.esa.int/api/objects','page[size]','100','page[number]','1',wo);
numpages = jsondata.meta.pagination.totalPages;

fprintf('Number of pages of data on discosweb for objects: %i \n', numpages);
fprintf('Downloading ~10 pages per minute \n');

clear d;
d = jsondata.data;

%% need to do 20 pages at a time, due to "TOO MANY REQUESTS" error 
for ind = 1:numpages
    if ind > 20 && mod(ind,20)==1
        pause(60);
        fprintf('\n');    
    end

    jsondata = webread('https://discosweb.esoc.esa.int/api/objects','page[size]','100','page[number]',num2str(ind),wo);
    d = [d; jsondata.data];
    fprintf('%i ',ind);    
end

fprintf('DISCOS data number of objects: %i \n', size(d,1))
fprintf('DISCOS data headers: \n')
disp(fieldnames(d(1).attributes))

% save('d_2023.mat','d')

% saved 8/27/2022: old/d_all.mat
% saved 3/10/2023: d_2023.mat

% DISCOS
% Treat as truth data for: NORAD ID, cosparID (launch yr), width, height, mass, object class

