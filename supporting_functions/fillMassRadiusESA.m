function [outmatsat,g1,g2,g3] = fillMassRadiusESA(inmatsat)
    % Data pipeline for mat_sats from initialized_01-2023.mat
    % 1) Get TLEs from space-track.org
    % 2) 
    % 3)
    % ESA method for filling in missing info from DISCOS
    % a) use space-track SATCAT designation for RCS size (S/M/L) and assign:
    %    0.1 / 1 / 10 m^2 area, for spherical density of aluminum (2,710 kg/m3)
    SATCATstname = 'spacetrack_satcat_03_2023.csv';
    satcatfn = which(SATCATstname);
    opts = detectImportOptions(SATCATstname); 
    satcatstdata = readtable(SATCATstname,opts);
    fprintf('Using RCS info from: %s \n', satcatfn);
    
    satcatSatnum = satcatstdata.NORAD_CAT_ID;

        % MATSATS DEFINITION
        idx_a = 1; idx_ecco = 2; idx_inclo = 3; 
        idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; 
        idx_bstar = 7; idx_mass = 8; idx_radius = 9;
        idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
        idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
        idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
        idx_objectclass = 23; idx_ID = 24;

    % objects without radius
    noRind = inmatsat(:,idx_radius) == 0;
    noRSatnum = inmatsat(noRind,idx_ID);
    
    [~,b,c] = intersect(noRSatnum,satcatSatnum);
    
    nlg = cellfun(@(x) strcmp(x,'LARGE'),satcatstdata(c,:).RCS_SIZE);
    nmd = cellfun(@(x) strcmp(x,'MEDIUM'),satcatstdata(c,:).RCS_SIZE);
    nsm = cellfun(@(x) strcmp(x,'SMALL'), satcatstdata(c,:).RCS_SIZE);

    fprintf('Objects without radius: %i small, %i medium, %i large objects \n ... now assigned with 0.1,1,10 m^2 area\n',...
        sum(nsm),sum(nmd),sum(nlg));
    fprintf('Summary: %i of %i objects without radius are now assigned \n',...
        sum(nsm + nmd + nlg), sum(noRind));
    
    outmatsat = inmatsat;
    noRinds = find(noRind);
    outmatsat(noRinds(nlg),idx_radius) = sqrt(10/pi);  % 10 m^2
    outmatsat(noRinds(nmd),idx_radius) = sqrt(1/pi);  % 1 m^2
    outmatsat(noRinds(nsm),idx_radius) = sqrt(0.1/pi);  % 0.1 m^2

    % categorize "unknown" objects that are still 0 radius    
%     noRSatnum(logical(nlg + nmd + nsm)) = [];  % find remainding noRsats
    noRind = outmatsat(:,idx_radius) == 0;
    noRSatnum = outmatsat(noRind,idx_ID);
    noRinds = find(noRind);

    [~,~,cc] = intersect(noRSatnum,satcatSatnum);

    % PAYLOAD -> LARGE
    % DEBRIS -> SMALL
    for ci = 1:numel(cc)
        c = cc(ci);
        curType = satcatstdata(c,:).OBJECT_TYPE{1};
        switch curType 
            case 'DEBRIS'
                outmatsat(noRinds(ci),idx_radius) = sqrt(0.1/pi);
            case 'PAYLOAD'
                outmatsat(noRinds(ci),idx_radius) = sqrt(10/pi);
            otherwise
                warning('SAT ID %i has type %s, skipping', satcatSatnum(c), curType);
        end
    end
    
    % ADD MASS (spherical aluminum for ESA; 2,710 kg/m3)
    noMind = outmatsat(:,idx_mass) == 0;
    outmatsat(noMind,idx_mass) = 4/3 * pi * outmatsat(noMind,idx_radius).^3 * 2710; % kg
    fprintf('Added mass to %i objects without mass\n', sum(noMind));

    [g1,g2,g3] = getZeroGroups(inmatsat); % fill out g1,g2,g3
    g1.gm = gmdistribution([sqrt(10/pi), 4/3 * pi * sqrt(10/pi).^3 * 2710],[0,0]); % Payload mu for [r,m]
    g2.gm = gmdistribution([sqrt(10/pi), 4/3 * pi * sqrt(10/pi).^3 * 2710],[0,0]); % RB 
    g3.gm = gmdistribution([sqrt(0.1/pi), 4/3 * pi * sqrt(0.1/pi).^3 * 2710],[0,0]); % Deb

end
