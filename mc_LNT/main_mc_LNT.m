%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT ARCLab's implementation of a space enivormental Simulation
%     MIT Orbital Capacity Tool - Monte Carlo (MOCAT-MC)
% 
% Authors: Richard Linares, Daniel Jang, Davide Gusmini, Andrea D'Ambrosio, 
%       Pablo Machuca, Peng Mun Siew
% https://github.mit.edu/arclab/orbitalrisk_MC
% 
% 08/23/2022: Initial commit
% 10/27/2022: Functionalized MOCAT-MC engine 
% 11/   2022: Vectorized subroutines
% 12/21/2022: main_mc.m created from main_mc_func_vec.m and main_mc_summary.m
% 12/03/2023: created from main branch main_mc_LNTiadc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_mc_LNT(MCconfig, RNGseed, LNTLC, Pfactor)
    % MCconfig is a structure of all relevant variables from a config setup file
    % usually in the 'configFiles' folder
  
    % Initialize RNG seed
    if (exist('RNGseed', 'var'))
        rng(RNGseed);
        fprintf('main_mc specified with seed %i \n', RNGseed);
    elseif isfield(MCconfig, 'seed')
        rng(MCconfig.seed);
        fprintf('main_mc specified with config seed %i \n', MCconfig.seed);
    end

    if ~exist('LNTLC', 'var')
        LNTLC = 0.1;  % 10 cm limit as default
    end
    fprintf('LNT characteristic length limit set to %0.3e meters \n', LNTLC);

    if ~exist('Pfactor', 'var')
        Pfactor = 1;  % no change to probability of collision calculation
    end
    fprintf('Pfactor set as %0.3e \n', Pfactor);

    addpath(genpath(pwd))

    % LOAD CFG 
    if ischar (MCconfig)
        MCconfig = eval(MCconfig);
    end
    loadCFG(MCconfig);  % see bottom    
    
    if isfield(MCconfig,'paramSSEM')
        param.paramSSEM = MCconfig.paramSSEM;        
        param.paramSSEM.dt_days = MCconfig.dt_days;
    end
    
    if isfield(MCconfig,'sample_params')
        param.sample_params = MCconfig.sample_params;
    else
        param.sample_params = 0;        
    end

    if ~exist('PMDconstel','var')
        PMDconstel = PMD;
        fprintf('PMDconstel set as PMD: %0.2e \n', PMDconstel);
    else
        fprintf('PMDconstel from cfg as %0.2e \n', PMDconstel);
    end

    getidx;

    % set dictionary for ConstellationID -> CompanyID
    k = [1;0];  v = [0;0];
    k = [k; unique(mat_sats(:,idx_constel)); unique(repeatLaunches(:,idx_constel))];
    v = [v; zeros(size(unique(mat_sats(:,idx_constel)))); zeros(size(unique(repeatLaunches(:,idx_constel)))) ];
    if ~isempty(constellation)
        k = [k; constellation.constIndex];
        v = [v; constellation.constCompanyID];
    end
    D = dictionary(k,v);

    
    % remove large data embedded in cfg (for saving)
    MCconfig.a_all = {};
    MCconfig.ap_all = {};
    MCconfig.aa_all = {};
    MCconfig.launchMC_step = {};

    % Assign constants to param structure for some functions
    param.req = radiusearthkm;
    param.mu = mu_const;
    param.j2 = j2;
    param.max_frag = max_frag;
    paramSSEM.species = [1,1,1,0,1,0];  % species: [S,D,N,Su,B,U];

    % Density profile
    param.density_profile = density_profile;
    
    % PREALLOCATE
    zmIdx = mat_sats(:,idx_mass) == 0;
    mat_sats(zmIdx,:) = [];
    zrIdx = mat_sats(:,idx_radius) == 0;
    mat_sats(zrIdx,:) = [];
    fprintf('mat_sats: %i zero-mass and %i zero-radius objects removed \n',...
        sum(zmIdx), sum(zrIdx));

    % setup repeatLaunches
    if (~exist('asemmat', 'var') || isempty(asemmat))
        dt = juliandate(time0) - min(repeatLaunches(:,idx_launch_date));    % align start of repeatLaunches to t0 [days]
        repeatLaunches(:,idx_launch_date) = repeatLaunches(:,idx_launch_date) + dt;
    end
    if launchRepeatSmooth
        % launchRepeatRangeDays = diff(launchRepeatYrs) * YEAR2DAY;
        launchRepeatRangeDays = julianFracYr(launchRepeatYrs(2)) - julianFracYr(launchRepeatYrs(1));
        repeatLaunches(:,idx_launch_date) = linspace( min(repeatLaunches(:,idx_launch_date)), ...
            min(repeatLaunches(:,idx_launch_date)) + launchRepeatRangeDays, size(repeatLaunches,1)); % linearly spaced
    end
    launchMC_step = repeatLaunches;  % rolling launch for 'matsat' and 'random_calendar'; initialize first cycle

    n_sats = size(mat_sats,1);
    fprintf('mat_sats: initial population: %i objects \n',n_sats);
    

    % Start diary if specified
    if (exist('save_diaryName', 'var') && ~isempty(save_diaryName))
        fidout = fopen([save_diaryName '.out'],'at');
        fiderr = fopen([save_diaryName '.err'],'at');
    else
        fidout = 0;
        fiderr = 0;
    end
        
    numObjects = zeros(n_time,1); % count the number of total object vs time
    numObjects(1,1) = n_sats;
    numObjTrigger = 1e6;  % save when this number of obj is reached & bump up as it's reached
    count_launch = uint16(zeros(n_time,1)); % count the number of collisions vs time
    count_coll = uint8(zeros(n_time,1)); % count the number of collisions vs time
    count_expl = uint8(zeros(n_time,1)); % count the number of explosions vs time
    count_debris_coll = uint8(zeros(n_time,1)); % number of collision debris vs time
    count_debris_expl = uint8(zeros(n_time,1)); % number of explosion debris vs time
    if save_output_file==3||save_output_file==4
        sats_info = cell(1,3);
    else
        sats_info = cell(n_time,3);         % contains info for SSEM binning
    end
    matsatsperN = {};
    frag_info = cell(n_time,4);         % contains info on cube statistics
    % frag_info5 = cell(n_time,1);        % contains info on collision statistics
    frag_info5vec  = [];                % collision info in vector
    frag_info6vec  = single([]);                % cube info in vector
    % count of species vs species in cube for all time
    frag_infovec = table('size',[n_time,16], 'VariableNames', ...
        {'SS','SC','SD','SN','SB','CC','CD','CN','CB','DD','DN','DB','NN','NB','BB','#above2percube'},...
        'VariableTypes',repelem({'uint32'},1,16));     % count of species vs species in cube for all time
    num_pmd = 0;
    num_deorbited = 0;
%     launch = 0;
%     out_future = [];
    count_tot_launches = 0;
    file_save_index = 0;
    S_MC = nan(n_time,numel(paramSSEM.R02)-1); D_MC = S_MC;  N_MC = S_MC; B_MC = S_MC;
    param_mean = nan((numel(paramSSEM.R02)-1)*4,3,n_time); param_var = param_mean;
    param_median = param_mean;
    
    % MATSATS DEFINITION
    % idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9;
    % idx_error = 10; idx_controlled = 11; idx_a_desired = 12; idx_missionlife = 13; idx_constel = 14; idx_date_created = 15; idx_launch_date = 16;
    % idx_r = [17 18 19]; idx_v = [20 21 22]; idx_objectclass = 23; idx_ID = 24;
    
    param.maxID = max([max(mat_sats(:,idx_ID)),0]);
    
    %launches_current_step_vec and launches_current_step_eq_vec
    idx_launch_in_extra = [idx_ID]; %additional columns to create for new launches
    
    %prop_mit_vec
    idx_prop_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_bstar,idx_controlled];
    idx_thalassa_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_bstar,idx_mass,idx_radius,idx_r,idx_v];
    %mat_sat_out = [a,ecco,inclo,nodeo,argpo,mo,errors,r_eci,v_eci];
    idx_prop_out = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_error,idx_r,idx_v];
    
    %orbcontrol_vec
    idx_control_in = [idx_a,idx_ecco,idx_inclo,idx_nodeo,idx_argpo,idx_mo,idx_controlled,idx_a_desired,idx_missionlife,idx_launch_date,idx_r,idx_v];
    %mat_sat_out = [a,controlled,r,v];
    idx_control_out = [idx_a,idx_controlled,idx_r,idx_v];
    
    %frag_exp_SBM_vec
    idx_exp_in = [idx_mass,idx_radius,idx_r,idx_v,idx_objectclass];
    
    %frag_col_SBM_vec
    idx_col_in = [idx_mass,idx_radius,idx_r,idx_v,idx_objectclass];

    objclassint_store = mat_sats(:,idx_objectclass);
    a_store = mat_sats(:,idx_a); 
    controlled_store = mat_sats(:,idx_controlled); 

    if save_output_file==3||save_output_file==4
        sats_info{1} = int8(objclassint_store);
        sats_info{2} = single(a_store);
        sats_info{3} = int8(controlled_store);    
        [popSSEM] = Fast_MC2SSEM_population(sats_info,paramSSEM);
        S_MC(1,:) = popSSEM(:,1);
        D_MC(1,:) = popSSEM(:,2);
        N_MC(1,:) = popSSEM(:,3);
        B_MC(1,:) = popSSEM(:,4);
    elseif save_output_file==6 
        [popSSEM,popSSEM_param_mean, popSSEM_param_var, popSSEM_param_median] = MC2SSEM_population_dist(mat_sats,paramSSEM);
        S_MC(1,:) = popSSEM(:,1);
        D_MC(1,:) = popSSEM(:,2);
        N_MC(1,:) = popSSEM(:,3);
        B_MC(1,:) = popSSEM(:,4);
        param_mean(:,:,1) = popSSEM_param_mean;
        param_var(:,:,1) = popSSEM_param_var;
        param_median(:,:,1) = popSSEM_param_median;
    elseif save_output_file >= 10
        % matsatsperN{1} = single(mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar,  idx_objectclass, idx_controlled]));
        % XXYYZ = p1_objectclass * 1000 + p1_const * 10 + p1_cont;
        XXYYZ = mat_sats(:,idx_objectclass) * 1000 + mat_sats(:,idx_constel) * 10 + mat_sats(:,idx_controlled);
        matsatsperN{1} = single([mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar, idx_ID]),XXYYZ]);
    else
        sats_info{1,1} = int8(objclassint_store);
        sats_info{1,2} = single(a_store);
        sats_info{1,3} = int8(controlled_store); 
    end
    
    [nS, nD, nN, nB] = categorizeObj(objclassint_store, controlled_store);
    if save_output_file >= 11
        savevars = {'MCconfig','param','paramSSEM','matsatsperN','frag_info5vec'};
        save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
    elseif save_output_file > 0    % test save (check for directory and write permission)
        save([filename_save(1:end-4),'_part_1.mat'],'-v7.3','param');
        fprintf('Test saving of %s successful\n', filename_save);
    end

    outmsg = sprintf('%i - %03i,\t PMD %04i,\t Deorbit %03i,\t Launches %03i,\t nFrag %03i,\t nCol %03i,\t nObj %i (%i,%i,%i,%i)\n', ...
        year(time0),day(time0,'dayofyear'), num_pmd ,num_deorbited, ...
        0, count_expl(1), count_coll(1), numObjects(1,1), nS, nD, nN, nB);
    fprintf(outmsg);
    if fidout > 0   % if diary is being written
        fprintf(fidout,outmsg);
    end


%START PROPAGATION
    for n = 2:n_time
        current_time = time0+days(tsince(n)/DAY2MIN);
        jd = juliandate(current_time);

    %PROPAGATION (one timestep at a time)
        n_sats = size(mat_sats,1);
        X_eci = zeros(n_sats,6);
        if ~exist('propagator','var')
            propagator = 'MIT';
        end
        
        if strcmpi(propagator,'SGP4')         % select between SGP4 or MIT's own propagator (deprecated at the moment)
            deorbit=[];
            for i = 1:n_sats
                [sats{i}, X_eci_temp, ~,~]=spg4_ecf(sats{i},tsince(n)); %tsince(n) tsince(n)

                % deorbit if spg4 outputs zero position vector or if satellite is
                % below 90 km altitude
                if norm(X_eci_temp(1,:))==0 || ((sats{i}.a*radiusearthkm)*(1-sats{i}.ecco)-radiusearthkm)<150 || norm(X_eci_temp(1:3,:))<(radiusearthkm+100) || sats{i}.error~=0 || sats{i}.a<0
                    deorbit=[deorbit;i];
                    X_eci(i,:)=0;    
                    sats{i}.r=[0 0 0];
                    sats{i}.v=[0 0 0];
                else
                    if norm(X_eci_temp(4:6,1))>30 % <<<<<<<<<<????
                        warning('satellite velocity > 30 km/s! ')
                    end
                    X_eci(i,:)=X_eci_temp;
                    sats{i}.r=X_eci_temp(1:3,1)';
                    sats{i}.v=X_eci_temp(4:6,1)';
                end
            end
        elseif strcmpi(propagator,'THALASSA') % use thalassa
            param.jd=jd;
            if n>1
                dt = 60*(tsince(n)-tsince(n-1));% units of time in seconds
            else
                dt = 60*tsince(n);% units of time in seconds
            end
            for i = 1:n_sats
                mat_sats(i,idx_prop_out) = prop_thalassa(mat_sats(i,idx_thalassa_in),dt,param);
            end
            deorbit = find(mat_sats(:,idx_r(1))==0 | (mat_sats(:,idx_a)*radiusearthkm).*(1-mat_sats(:,idx_ecco))<(150+radiusearthkm) | sqrt(sum(mat_sats(:,idx_r).^2,2))<(radiusearthkm+100) | mat_sats(:,idx_error)~=0 | mat_sats(:,idx_a)<0);           
        else % use prop_mit 
            param.jd=jd;
            
            %prop_mit_vec
            if n>1
                [mat_sats(:,idx_prop_out)] = prop_mit_vec(mat_sats(:,idx_prop_in),60*(tsince(n)-tsince(n-1)),param);% units of time in seconds
            else
                [mat_sats(:,idx_prop_out)] = prop_mit_vec(mat_sats(:,idx_prop_in),60*tsince(n),param);% units of time in seconds
            end

            % deorbit if spg4 outputs zero position vector or if satellite is
            % below 90 km altitude
            deorbit = find(mat_sats(:,idx_r(1))==0 | (mat_sats(:,idx_a)*radiusearthkm).*(1-mat_sats(:,idx_ecco))<(150+radiusearthkm) | sqrt(sum(mat_sats(:,idx_r).^2,2))<(radiusearthkm+100) | mat_sats(:,idx_error)~=0 | mat_sats(:,idx_a)<0);           
        end     
        % Remove all deorbited objects
        num_deorbited = length(deorbit);
        mat_sats(deorbit,:) = [];
        
    %ORBIT CONTROL        
        % if mod(n,step_control)==1 || step_control==1
            %orbcontrol_vec
            % [mat_sats(:,idx_control_out),deorbit_PMD] = orbcontrol_vec(mat_sats(:,idx_control_in),tsince(n),time0,orbtol,PMD,DAY2MIN,YEAR2DAY,param);
            [mat_sats(:,idx_control_out),deorbit_PMD, deorbitfail] = ...
                orbcontrol_constel(mat_sats(:,[idx_control_in,idx_constel,idx_objectclass]),...
                tsince(n),time0,orbtol,PMD,PMDconstel,DAY2MIN,YEAR2DAY,param);
            
             % RB PMD if cfg.PMDrb exists
             if exist('PMDrb','var')
                 find_RBs = find(mat_sats(:,idx_objectclass) == 5);
                 find_RBrecent = find_RBs(mat_sats(find_RBs,idx_launch_date) > juliandate(current_time) - dt_days);
                 rand_life = rand(numel(find_RBrecent),1);
                 check_PMD = PMDrb < rand_life; %check if PMD is fulfilled (fails)
                 deorbitRB = find_RBrecent(~check_PMD); % successful PMDs: deorbit
                 deorbit_PMD = [deorbit_PMD; deorbitRB];
                 if ~isempty(check_PMD)
                     fprintf('\t\t %i/%i RB PMD \n', sum(~check_PMD), numel(check_PMD));
                 end
             end
             % save constellation ID for deorbited sats for launches later
             deorbitConstels = mat_sats([deorbit_PMD;deorbitfail], idx_constel);
            
            % Remove all post-mission disposal satellites
            num_pmd = length(deorbit_PMD);
                % DEBUG: AGE and SMA OF DEORBITING SATS
%                 figure(10); subplot(211); histogram((mat_sats(deorbit_PMD,idx_launch_date) - jd)/YEAR2DAY); xlim([-10,0]); xlabel('years old'); title('PMD')
%                 subplot(212); histogram((mat_sats(deorbit_PMD,idx_a) -1) * radiusearthkm, 200:50:2000); xlim([200, 2000]); xlabel('Altitude (km)')
            mat_sats(deorbit_PMD,:) = [];
        % else
        %     num_pmd = 0;
        % end
    
    %EXPLOSIONS (for Rocket Body OR PAYLOADS)
        n_sats = size(mat_sats,1);
        out_frag = [];
        % find_rocket = find(mat_sats(:,idx_objectclass)==5); %Rocket bodies
        find_rocket = find(mat_sats(:,idx_objectclass)==5 | mat_sats(:,idx_objectclass)==1); %Rocket bodies OR PAYLOADS
        rand_P_exp = rand(numel(find_rocket),1); %random numbers for each rocket body
        if exist('P_frag_cutoff','var')
            find_P_exp = find(rand_P_exp<P_frag & ...
                (jd - jd2date(mat_sats(find_rocket,idx_launch_date))) < P_frag_cutoff);
        else
            find_P_exp = find(rand_P_exp<P_frag);
        end
        remove_frag = find_rocket(find_P_exp); %Pre-identify objects to remove
        for idx_P_exp_temp = numel(find_P_exp):-1:1 %reverse order in case elements need to be deleted from remove_frag variable
            idx_P_exp = remove_frag(idx_P_exp_temp);
            
            p1_all = mat_sats(idx_P_exp,:);
            p1_mass = p1_all(1,idx_mass);
            p1_objectclass = p1_all(1,idx_objectclass);
            p1_constel = p1_all(1,idx_constel);
            p1_in = p1_all(1,idx_exp_in);

            %frag_exp_SBM_vec
            %p1_in = [p1.mass,p1.radius,p1.r,p1.v,p1.objectclass]
            [debris1] = frag_exp_SBM_vec(tsince(n), p1_in, param);
            param.maxID = param.maxID+size(debris1,1); %update maximum ID in population, since new objects have been added
                        
            if isempty(debris1) %if no debris is generated, do not remove object from simulation
                remove_frag(idx_P_exp_temp) = [];
            else %if some debris is generated, add it to out_frag matrix
                debris1(:,idx_constel) = p1_constel;  % add parent constellation (forensics)
                out_frag = [out_frag; debris1];
                XXYYZ1 = p1_objectclass * 1000 + p1_const * 10 + p1_cont;
                errmsg = sprintf('Year %i - Day %03i \t Explosion, p1 XXYYZ %i, %0.1f kg, nDebris %i\n', ...
                    year(current_time), day(current_time,'dayofyear'), XXYYZ1, p1_mass, size(debris1,1));
                fprintf(2,errmsg);
                if fiderr > 0  % DIARY
                    fprintf(fiderr,errmsg);
                end
                count_expl(n) = count_expl(n)+1;      
            end
        end
        mat_sats(remove_frag,:) = [];

    %COLLISIONS
        if ( exist('skipCollisions','var') && skipCollisions == 1 ) || isempty(mat_sats)
            collision_array = [];
        else
            collision_cell = cube_vec_v3(mat_sats(:,idx_r), CUBE_RES, collision_alt_limit);          % do CUBE
            collision_array = cell2mat(collision_cell);
        end 

        remove_collision = [];
        out_collision = [];
        P = [];
        if ~isempty(collision_array)
            p1_idx = collision_array(:,1);
            p2_idx = collision_array(:,2);
            p1_all = mat_sats(p1_idx,:);
            p2_all = mat_sats(p2_idx,:);


            p1_controlled = p1_all(:,idx_controlled); p1_radius = p1_all(:,idx_radius); p1_v = p1_all(:,idx_v); p1_const = p1_all(:,idx_constel); p1_comp = D(p1_const);
            p2_controlled = p2_all(:,idx_controlled); p2_radius = p2_all(:,idx_radius); p2_v = p2_all(:,idx_v); p2_const = p2_all(:,idx_constel); p2_comp = D(p2_const); 
            % probability of collision
            Pij = collision_prob_vec(p1_radius, p1_v, p2_radius, p2_v, CUBE_RES) ./ Pfactor;
                % probability of collision over dt: 5 days
                % P = 4/3*pi*((p1.radius^3*p2.radius^3)/CUBE_RES^3)*Pij*dt_days*DAY2SEC;
            P = zeros(size(p1_controlled));
            sum_controlled = p1_controlled + p2_controlled;
            check_0 = sum_controlled<0.5;

        % NEWLY ADDED for CONSTELLATIONS 
        f1 = sum_controlled > 0.5 & sum_controlled < 1.5; % just one controlled
        f2 = sum_controlled > 1.5;                        % both controlled
        f0constel = (p1_const + p2_const) == 0;           % neither in constel
        % simple "same constel" definition:
            % fsameconstel = (p1_const == p2_const) & (p1_const > 1); % strictly same constel (note: constel=1 means misc constel)
            % fdiffconstel = ((p1_const > 0) & (p2_const > 0) & (p1_const ~= p2_const)) | ((p1_const == 1) & (p2_const == 1));
        % same company definition:
            fsameconstel = (p1_comp == p2_comp) & (p1_comp > 0); % same company means same constel (comp 0 means no company and no constel)
            fdiffconstel = ((p1_comp > 0) & (p2_comp > 0) & (p1_comp ~= p2_comp)) | ((p1_const == 1) & (p2_const == 1)); % constellation 1 vs 1 = different constellation

            P(f1) = Pij(f1)*(alph*dt_days*DAY2SEC); % one controlled: alph
            if exist('alphlnt','var')
                f1lntidx = f1 & (p1_radius < 0.05 | p2_radius < 0.05); % 5 cm radius
                P(f1lntidx ) = Pij(f1lntidx )*(alphlnt*dt_days*DAY2SEC); % one controlled, but against LNT: alphlnt
            end
            P(f2 & f0constel) = Pij(f2 & f0constel)*(alph_a*dt_days*DAY2SEC);  % both controlled, neither constel: alph_a
            P(f2 & fsameconstel) = Pij(f2 & fsameconstel)*(alph_c_inter*dt_days*DAY2SEC); % same constel: alph_c_inter
            P(f2 & fdiffconstel) = Pij(f2 & fdiffconstel)*(alph_c_intra*dt_days*DAY2SEC); % diff constel: alph_c_intra
            P(check_0) = Pij(check_0)*(dt_days* DAY2SEC); % neither controlled
        
            % Monte-carlo simulation - there is a collision if random number is lower than prob of collision
            rand_P = rand(size(p1_controlled));
            find_P = find(rand_P<P);
            frag4tmp = zeros(numel(find_P),3);    % col4:   index for collision, debAn, debBn               [Mx3] uint16

            for idx_P_temp = 1:numel(find_P)
                idx_P = find_P(idx_P_temp);
                curP = P(idx_P);                
                p1_mass = p1_all(idx_P,idx_mass);
                p2_mass = p2_all(idx_P,idx_mass);
                p1_objectclass = p1_all(idx_P,idx_objectclass);
                p2_objectclass = p2_all(idx_P,idx_objectclass);
                p1_radius = p1_all(idx_P,idx_radius);
                p2_radius = p2_all(idx_P,idx_radius);
                p1_cont = p1_all(idx_P,idx_controlled);
                p2_cont = p2_all(idx_P,idx_controlled);
                p1_const = p1_all(idx_P,idx_constel);
                p2_const = p2_all(idx_P,idx_constel);
                p1_r = p1_all(idx_P,idx_r);  p2_r = p2_all(idx_P,idx_r);
    
                p1_in = p1_all(idx_P,idx_col_in); 
                p2_in = p2_all(idx_P,idx_col_in); 
    
                [debris1, debris2, iscata] = frag_col_SBM_vec_lc(tsince(n), p1_in, p2_in, param, LNTLC);

                param.maxID = param.maxID+size(debris1,1)+size(debris2,1); %update maximum ID in population, since new objects have been added
    
                out_collision = [out_collision; debris1; debris2];          % add debris
    
                if ~isempty(debris1) || ~isempty(debris2)
                    XXYYZ1 = p1_objectclass * 1000 + p1_const * 10 + p1_cont;
                    XXYYZ2 = p2_objectclass * 1000 + p2_const * 10 + p2_cont;
                    errmsg = sprintf('Year %i - Day %03i \t Collision, p1 XXYYZ %i, %0.1f kg, p2 %i, %0.1f kg, nDebris1 %i, nDebris2 %i, altkm %0.0f, curP %0.2e, isCata: %i\n', ...
                        year(current_time), day(current_time,'dayofyear'), ...
                        XXYYZ1, p1_mass, XXYYZ2,  p2_mass, size(debris1,1), size(debris2,1), ...
                        mean([norm(p1_r), norm(p2_r)])-radiusearthkm , curP, iscata);
                    fprintf(2,errmsg);
                    if fiderr > 0
                        fprintf(fiderr,errmsg);
                    end
                    count_coll(n) = count_coll(n)+1;
                    remove_collision = [remove_collision; p1_idx(idx_P); p2_idx(idx_P)]; % remove parent objects    
                end
                % col4:   index for collision, debAn, debBn               [Mx3] uint16
                frag4tmp(idx_P_temp,1) = idx_P;
                frag4tmp(idx_P_temp,2) = size(debris1,1);
                frag4tmp(idx_P_temp,3) = size(debris2,1);
                frag5tmpvec = [n, iscata, p1_mass, p2_mass, (mean([norm(p1_r), norm(p2_r)])-radiusearthkm), size(debris1,1) + size(debris2,1), p1_radius, p2_radius, p1_objectclass, p2_objectclass, p1_cont, p2_cont];
                frag_info5vec = [frag_info5vec; frag5tmpvec];
            end
        end


    %NEW LAUNCHES
    out_future = [];
        if strcmpi(launch_model,'matsat')  % for TLE.. or otherwise (mat_sats repeater)
        % elseif strcmpi(launch_model,'matsat')
            [out_future,launch,launchMC_step] = launches_current_step_eq_vec(launch_model,time0,current_time,tsince,n,DAY2MIN,YEAR2DAY,param,idx_launch_in_extra,...
                launchMC_step , repeatLaunches, launchRepeatYrs);
            % repeatLaunches: static list from cfg to repeat over -- currently setup for yearly repeat
        elseif strcmpi(launch_model,'random') 
            [out_future,launch] = launches_current_step_vec(launch_model,time0,current_time,tsince,n,param,idx_launch_in_extra,...
                total_launch_per_year * (1+launch_increase_per_year)^floor(years(current_time - time0)), dt_days);
        elseif strcmpi(launch_model,'data') || contains(launch_model,'Somma')
            [out_future,launch] = launches_current_step_eq_vec(launch_model,time0,current_time,tsince,n,param,idx_launch_in_extra,...
                launchMC_step,ind_launch,additional_launches,ind_launch_add);
        elseif strcmpi(launch_model,'random_calendar') 
            if n<n_time
                [out_future,launch,launchMC_step] = launches_current_step_eq_vec(launch_model,time0,current_time,tsince,n,DAY2MIN,YEAR2DAY,param,idx_launch_in_extra,...
                    launchMC_step);
            else
                out_future = [];
            end
        end

        % CONSTELLATION launches
        for ind = 1:size(constellation,1)
            jdFirstLaunch = julianFracYr(constellation.FirstLaunch(ind));
            jdEndLaunch = julianFracYr(constellation.FinishLaunch(ind));
            if isnan(constellation.EndOperation(ind))
                constellation.EndOperation(ind) = 5000; % far future [yr]
            end
            jdEndOperation = julianFracYr(constellation.EndOperation(ind));
            curConstID = constellation.constIndex(ind);
            curNramp = (constellation.TotalSatsPlanned(ind) - constellation.SatsOnStn(ind)) ...
                    / (jdEndLaunch - jdFirstLaunch) * YEAR2DAY;
            if (jd >= jdFirstLaunch) && (jd < jdEndOperation) && (jd < jdEndLaunch)
                % if still in ramp-up phase
                Nrampperyear(ind) = curNramp;
                Nmaintperyear(ind) = floor( (jd - jdFirstLaunch) / YEAR2DAY / constellation.missionlife(ind) ) * ...
                                        Nrampperyear(ind);  
            elseif  jd >= jdFirstLaunch && jd < jdEndOperation
                % if ramping is finished (constellation filled up), wait until deorbits are happening,
                % and replenish as they deorbit (maintenance mode)
                Nrampperyear(ind) = 0;
                if jd > jdFirstLaunch + YEAR2DAY * constellation.missionlife(ind)
                    % Nmaintperyear(ind) = constellation.TotalSatsPlanned(ind) / constellation.missionlife(ind);
                    Nmaintperyear(ind) = sum(deorbitConstels == curConstID) / dt_days * YEAR2DAY;
                else % still waiting for deorbits to happen
                    Nmaintperyear(ind) = 0;
                end

                % OR just count constIDs from missionLife ago
                
            else % before first launch, or after end of operation
                Nrampperyear(ind) = 0;
                Nmaintperyear(ind) = 0;
            end
        end
        % fprintf('ramp per year: %0.1f \t maint per year: %0.1f \n', Nrampperyear, Nmaintperyear)


        if size(constellation,1) > 0
            % debug
            fprintf('\t\t Num constellations Ramping up: %i \t\t Maintaining: %i \t\t Yet to start: %i \t\t Ended operation: %i \n', ...
                sum( jd >= julianFracYr(constellation.FirstLaunch) & jd <= julianFracYr(constellation.FinishLaunch)), ...
                sum( jd > julianFracYr(constellation.FinishLaunch) & jd <= julianFracYr(constellation.EndOperation)), ...
                sum( jd < julianFracYr(constellation.FirstLaunch)), sum( jd > julianFracYr(constellation.EndOperation)) );
            avglaunch = (Nmaintperyear + Nrampperyear)/(YEAR2DAY/dt_days);  % avg launch per timestep (fractional)
            intlaunch = floor(avglaunch);  % guaranteed launch amount per time step (int)
            addlaunch = rand(size(avglaunch)) < (avglaunch - intlaunch);    % additional launches - sampled from fractional part
            curConstLaunch = intlaunch + addlaunch;

            constMatsat = nan(size(avglaunch,2), 24);
            constMatsat(:,[idx_a, idx_ecco, idx_inclo, idx_constel, idx_bstar, idx_mass, idx_radius, idx_missionlife]) = ...
                [constellation.Altitude / radiusearthkm + 1, ones(size(avglaunch))'/1e6, ...
                deg2rad(constellation.Inclination), constellation.constIndex, ...
                0.5 * 2.2 * constellation.radius .^2 ./ constellation.mass * 0.157, ...
                constellation.mass, constellation.radius, constellation.missionlife];
            constMatsat(:, idx_a_desired) = constMatsat(:, idx_a);
            constMatsat(:,[idx_controlled, idx_objectclass]) = [ ones(size(avglaunch))',  ones(size(avglaunch))'];

            constLaunch = constMatsat(repelem(1:size(constMatsat,1), curConstLaunch), :  ); % populate constellation launch matsats

            % Scramble launch time within this time step (dt_days), argpo, mo, nodeo
            constLaunch(:,idx_launch_date) = jd + rand(size(constLaunch,1),1) * dt_days;
            constLaunch(:,[idx_argpo, idx_mo, idx_nodeo]) = 2 * pi * rand(size(constLaunch,1),3);

            param.maxID = param.maxID+size(out_future,1)+size(constLaunch,1);       %update maximum ID in population
            count_tot_launches = count_tot_launches+size(out_future,1)+size(constLaunch,1);
            out_future = [out_future; constLaunch];
        end

        count_launch(n) = size(out_future,1);
    

    %DATA PROCESSING
        % note this is put at the end because we assume these object don't
        % generate collisions during their first timestep of existance
        mat_sats(remove_collision,:) = [];                
        mat_sats = [mat_sats; out_future; out_frag; out_collision];
                % DEBUG: AGE and SMA OF DEORBITING SATS
%                 vizMatsats(out_future,1); vizMatsats(mat_sats,2); subplot(311); ylim([0,30]); 

    %ACCOUNTING
        n_sats = size(mat_sats,1);
        numObjects(n,1) = n_sats;
        
        objclassint_store = mat_sats(:,idx_objectclass);
        a_store = mat_sats(:,idx_a); 
        controlled_store = mat_sats(:,idx_controlled); 
        
        if save_output_file==3||save_output_file==4
            sats_info{1} = int8(objclassint_store);
            sats_info{2} = single(a_store);
            sats_info{3} = int8(controlled_store);    
            [popSSEM] = Fast_MC2SSEM_population(sats_info,paramSSEM); % [s_count' d_count' n_count' su_count' b_count' u_count' ];
            S_MC(n,:) = popSSEM(:,1);
            D_MC(n,:) = popSSEM(:,2);
            N_MC(n,:) = popSSEM(:,3);
            B_MC(n,:) = popSSEM(:,4);
        elseif save_output_file==6 
            [popSSEM,popSSEM_param_mean, popSSEM_param_var, popSSEM_param_median] = MC2SSEM_population_dist(mat_sats,paramSSEM);
            S_MC(n,:) = popSSEM(:,1);
            D_MC(n,:) = popSSEM(:,2);
            N_MC(n,:) = popSSEM(:,3);
            B_MC(n,:) = popSSEM(:,4);
            param_mean(:,:,n) = popSSEM_param_mean;
            param_var(:,:,n) = popSSEM_param_var;
            param_median(:,:,n) = popSSEM_param_median;
        end
        if save_output_file >= 10 % | save_output_file == 11
            if mod(n,saveMSnTimesteps) == 1
                % matsatsperN{1 + floor(n/saveMSnTimesteps)} = single(mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar,  idx_objectclass, idx_controlled]));
                XXYYZ = mat_sats(:,idx_objectclass) * 1000 + mat_sats(:,idx_constel) * 10 + mat_sats(:,idx_controlled);
                matsatsperN{1 + floor(n/saveMSnTimesteps)} = single([mat_sats(:,[idx_a, idx_ecco, idx_mass, idx_radius, idx_bstar, idx_ID]),XXYYZ]);
            end
        else
            sats_info{n,1} = int8(objclassint_store);
            sats_info{n,2} = single(a_store);
            sats_info{n,3} = int8(controlled_store);
        end
       
        if ~isempty(collision_array) && (save_output_file == 3 | save_output_file == 11)
                % COUNT COLLISION BY OBJECT TYPES
                cubeclasses = [p1_all(:,idx_objectclass), p2_all(:,idx_objectclass)];
                cubeconstel = [p1_all(:,idx_constel),   p2_all(:,idx_constel)];
                cubecontrol = [p1_all(:,idx_controlled), p2_all(:,idx_controlled)];
                % change obj types: S:1, Const:2, D:3, N:4, B:5
                fraginfoclass = zeros(size(cubeclasses));
                fraginfoclass(cubeclasses == 1 & cubeconstel <1 & cubecontrol == 1) = 1; % S
                fraginfoclass(cubeconstel >0) = 2;                     % Constel
                fraginfoclass(cubeclasses == 1 & cubecontrol ~= 1) = 3; % D
                fraginfoclass(cubeclasses == 5) = 5;                    % B
                fraginfoclass(cubeconstel <1) = 4;                    % all others: N
                fraginfoclass = sort(fraginfoclass,2);
                fraginfoclasspair = fraginfoclass(:,1) * 10 + fraginfoclass(:,2);
                colpaircount = [sum(fraginfoclasspair == 11), sum(fraginfoclasspair == 12),...
                    sum(fraginfoclasspair == 13), sum(fraginfoclasspair == 14),...
                    sum(fraginfoclasspair == 15), sum(fraginfoclasspair == 22),...
                    sum(fraginfoclasspair == 23), sum(fraginfoclasspair == 24),...
                    sum(fraginfoclasspair == 25), sum(fraginfoclasspair == 33),...
                    sum(fraginfoclasspair == 34), sum(fraginfoclasspair == 35),...
                    sum(fraginfoclasspair == 44), sum(fraginfoclasspair == 45),...
                    sum(fraginfoclasspair == 55), size(collision_array,1) - size(collision_cell,1)];
                       % SS2, SC3, SD4, SN5, SB6, 
                       % CC4, CD, CN, CB
                       % DD, DN, DB, NN, NB, BB, above 2 per cube
            frag_infovec{n,:} = colpaircount; 
            

        elseif ~isempty(collision_array) && (save_output_file == 4 || save_output_file == 5 || save_output_file == 7)
            % save fragmentation info:  
            % col1: P, p1 satno, p2 satno                           [Nx3] single
            % col2: p1 mass, p2 mass (0.1 kg)                       [Nx2] uint16
            % col3: objclass, objclass, dv (0.1 km/s), alt10km      [Nx4] uint8
            % col4: index for collision, debAn, debBn               [Mx3] uint16
            % "col5": p1 mass, p2 mass (0.1 kg), alt10km, debAn, debBn [Mx5] uint16
            frag_info{n,1} = single([P, p1_all(:, idx_ID), p2_all(:, idx_ID)]);
            frag_info{n,2} = uint16([p1_all(:, idx_mass), p2_all(:, idx_mass)] * 10);
            frag_info{n,3} = uint8([p1_all(:, idx_objectclass), p2_all(:, idx_objectclass), ...
                                vecnorm((p1_all(:,idx_v) - p2_all(:,idx_v))')' * 10, ...
                                (mean([vecnorm(p1_all(:,idx_r),2,2), vecnorm(p2_all(:,idx_r),2,2)],2)-radiusearthkm)/10 ...
                                ]);
            frag_info{n,4} = single(frag4tmp);
%             frag_info5{n} = uint16(frag5tmp);
        end
            
    
        count_debris_coll(n,1) = size(out_collision,1);
        count_debris_expl(n,1) = size(out_frag,1);

        [nS, nD, nN, nB] = categorizeObj(objclassint_store, int8(controlled_store));
        outmsg = sprintf('%i-%03i,\t PMD %02i,\t Deorbit %03i,\t Launches %03i/%05i,\t nExp %02i/%04i,\t nCol %02i/%04i,\t EV(col) %0.2e, \t nObj %i (%i,%i,%i,%i)\n', ...
            year(current_time), day(current_time,'dayofyear'), num_pmd ,num_deorbited, ...
            count_launch(n),sum(count_launch), count_expl(n), sum(count_expl(1:n)), ...
            count_coll(n), sum(count_coll(1:n)), sum(P), numObjects(n,1), nS, nD, nN, nB);
        % fprintf(outmsg);
        fprintf('popzero: %i \t %s',sum( mat_sats(:,idx_ID) == -1), outmsg);
        if fidout > 0
            fprintf(fidout, outmsg);
        end

    % CHECK MEMORY CONSUMPTION
        sat_bytes = whos('sats_info').bytes;
        frag_info_bytes = whos('frag_info').bytes;
        % frag_info5_bytes = whos('frag_info5').bytes;
        frag_info6_bytes = whos('frag_info6vec').bytes;
        matsatsperN_bytes = whos('matsatsperN').bytes;
        total_bytes = sat_bytes + frag_info_bytes + frag_info6_bytes + matsatsperN_bytes; % + frag_info5_bytes 

    % SAVE OUTPUT 
        if save_output_file > 0 && ( (mod(n,n_save_checkpoint) == 1 || n==n_time || numObjects(n,1) > numObjTrigger || total_bytes/(1e+9) > 30) )
            if total_bytes/(1e+9) > 30
                warning('Total file size exceeded 30GB. Total file size is %s GB.',total_bytes/(1e+9));
            end

            if numObjects(n,1) > numObjTrigger
                numObjTrigger = numObjTrigger + 5e5;  % increase limit
            end
            
            % Increment file save index by 1
            file_save_index = file_save_index +1;

            switch save_output_file
                case 1  % save entire workspace
                    savevars = {'*'}; 
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 2  % save sats_info, config, params
                    savevars = {'sats_info','MCconfig','param','paramSSEM'};
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 3  % save summary (S_MC, N_MC, D_MC, etc)
                    savevars = {'S_MC','D_MC','N_MC','B_MC','MCconfig','param','paramSSEM','frag_info5vec','frag_infovec'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 4  % save summary and collision stats
                    savevars = {'S_MC','D_MC','N_MC','B_MC','MCconfig','param','paramSSEM','frag_info','frag_info5vec','frag_infovec'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 5  % save just collision stats
                    savevars = {'MCconfig','param','paramSSEM','frag_info','frag_info5vec','frag_infovec'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 6  % Tracking mean/var/median of physical parameters
                    savevars = {'S_MC','D_MC','N_MC','B_MC','MCconfig','param','paramSSEM','param_mean','param_var','param_median'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 7  % Collision tracking (UROP)        
                    frag_info7 = cell(n_time,2); % just 2 columns needed
                    frag_info7(:,1) = frag_info(:,1);  % 1st col
                    frag_info7(:,2) = frag_info(:,4);  % 4th col
                    savevars = {'frag_info7'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                    clear frag_info7;
                case 10  % save all, decimate matsats
                    savevars = {'MCconfig','param','paramSSEM','matsatsperN'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                case 11  % save decimated matsats & all fraginfo5
                    savevars = {'MCconfig','param','paramSSEM','matsatsperN','frag_info5vec','frag_infovec','frag_info6vec'};
                    MCconfig.mat_sats = [];  % to save space
                    save([filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat'],'-v7.3',savevars{:});
                otherwise
                    error('save_output_file flag set to unsupported value: %i', save_output_file);
            end
    
            if save_output_file==3||save_output_file==4
                sats_info = cell(1,3);
            else
                sats_info = cell(n_time,3);         % contains info for SSEM binning
            end
            frag_info6vec  = single([]);        % cube pair info in vector
            frag_info = cell(n_time,4);         % contains info on cube statistics
            % frag_info5 = cell(n_time,1);        % contains simplified info on collision statistics
            warning('Large variables cleared and %s saved', [filename_save(1:end-4),'_part_',num2str(file_save_index),'.mat']);
        end
    end
    
    fclose all;
    sprintf('\n === FINISHED MC RUN (main_mc.m) WITH SEED: %i === \n');

end

function loadCFG(cfg)
    fns = fieldnames(cfg);
    for ind = 1:numel(fns)
        if strcmp(fns{ind}, 'param')
            param = cfg.param;
            assignin('caller','param', param);
        else
            assignin('caller',fns{ind}, cfg.(fns{ind}));
        end
    end
end

function [nS, nD, nN, nB] = categorizeObj(objint_st, cont_st)
    % definitions from Fast_MC2SSEM_population.m
    nS = sum(objint_st==1 & cont_st==1);
    nD = sum(objint_st==1 & cont_st~=1);
    nN = sum(objint_st==3 | objint_st==4 | objint_st>=7);
    nB = sum(objint_st==5 | objint_st==6);
end