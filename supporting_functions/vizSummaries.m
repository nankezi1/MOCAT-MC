function out = vizSummaries(varargin)
% Plots S/D/N/B/total populations, collision stats of all relevant summary files
% input can be folders (./0728_ADEPT/) wildcards (./0728_ADEPT/su*.mat)

% frag5tmpvec = [n, iscata, p1_mass, p2_mass, (mean([norm(p1_r), norm(p2_r)])-radiusearthkm), size(debris1,1) + size(debris2,1), p1_radius, p2_radius, p1_objectclass, p2_objectclass, p1_cont, p2_cont];

% close(10:12);

fs = [];
for n = 1:nargin
    curarg = varargin{n};
    if endsWith(curarg, '.mat')
        f = dir(curarg);
        fs = [fs; f];
    end
end

disp({fs.name})

for n = 1:numel(fs)
    [~, parentfolder, ~] = fileparts(fs(n).folder);
    fprintf('Folder: %s \t File: %s \n', parentfolder, fs(n).name);

    clear('MCconfig','S_MC','N_MC','B_MC','D_MC','fi5summary','fi5summ','paramSSEM')
    load(fullfile(fs(n).folder, fs(n).name),'MCconfig','S_MC','N_MC','B_MC','D_MC','fi5summary','fi5summ','paramSSEM');
    
    if exist('fi5summary','var')
        fi5summary(cellfun(@isempty,(fi5summary)))  = [];  % clear empty fi5summary
    elseif exist('fi5summ','var')
        % fi5summary  = cell2mat(fi5summ');
        fi5summary  = fi5summ;
    else
        fi5summary = {};
    end

    % vizSDNs(S_MC,D_MC,N_MC,B_MC,paramSSEM,MCconfig);

    if isfield(MCconfig,'saveMSnTimesteps')
        ts = [1:size(S_MC,1)] * MCconfig.saveMSnTimesteps * MCconfig.dt_days / MCconfig.YEAR2DAY;  % yrs
    else
        ts = [1:size(S_MC,1)] * MCconfig.dt_days / MCconfig.YEAR2DAY;  % yrs
    end

    % total number
    figure(100); hold on;
        totnum = squeeze(sum(S_MC,2)) + squeeze(sum(D_MC,2)) + ...
             squeeze(sum(N_MC,[2,3]));
        if ~isempty(B_MC)
            totnum = totnum + squeeze(sum(B_MC,2));
        end
        p = plot(ts,mean(totnum,2),'LineWidth',2,'DisplayName',[parentfolder filesep fs(n).name]);
        hold on;
        plot(ts,totnum','Color',[p.Color, 0.2],'HandleVisibility','off');

        title('Total number of objects'); 
        ylabel('Number'); xlabel('Years');
        xlim([0,ts(end)]);

    % Number of Payloads
    figure(110); 
        plot(ts, mean(sum(S_MC,2),3),'DisplayName',[parentfolder filesep fs(n).name]);
        hold on;
        title('Payload count');
        ylabel('Number'); xlabel('Years');
        xlim([0,ts(end)]);
        a = gca; ylim([0,a.YLim(2)])


    % total number vs altitude (imagesc)
        tot_MC = S_MC + D_MC + B_MC + squeeze(sum(N_MC,3));
        tot_MC_mean = mean(tot_MC,3);
        nanIdx = isnan(tot_MC_mean(:,1)); % in time
        valid_tot_mean = tot_MC_mean(~nanIdx,:); % w/o nan'd times

    figure(1000+n); imagesc(ts(~nanIdx), paramSSEM.R02, valid_tot_mean');
        a = gca; a.YDir = 'normal'; xlabel('Years'); ylabel('Altitude (km)');
        title('Total Population',fs(n).name,'interpreter','none'); colormap('jet'); 
        colorbar;

    % cumulative collisions
    if ~isempty(fi5summary)
    figure(101); hold on;
        % create time boundary bins to count collisions (yrs)
        totcol = cellfun(@(x) histcounts(x(:,1) * MCconfig.dt_days / MCconfig.YEAR2DAY , ...
            [0, ts] )', fi5summary, 'UniformOutput', false);
        totcol = cell2mat(totcol);  % nT x nMC

        p = plot(ts,mean(cumsum(totcol),2),'LineWidth',2,'DisplayName',[parentfolder filesep fs(n).name]);
        hold on;
        plot(ts,cumsum(totcol)','Color',[p.Color, 0.2],'HandleVisibility','off');

        title('Cumulative number of collisions');
        ylabel('Number'); xlabel('Years')
        xlim([0,ts(end)]); grid on;
    

    % collisions per year
    figure(102); hold on;
        % create time boundary bins to count collisions (yrs)
        totcol = cellfun(@(x) histcounts(x(:,1) * MCconfig.dt_days / MCconfig.YEAR2DAY , ...
            [0, ts] )', fi5summary, 'UniformOutput', false);
        totcol = cell2mat(totcol);  % nT x nMC

        p = plot(ts,mean(totcol,2),'LineWidth',2,'DisplayName',[parentfolder filesep fs(n).name]);
        hold on;
        plot(ts,totcol','Color',[p.Color, 0.2],'HandleVisibility','off');

        title('Collisions per year');
        ylabel('Number'); xlabel('Years')
        xlim([0,ts(end)]); grid on;

    % Collision per altitude (across all duration) 
    figure(103); hold on;
        altbins = 150:50:2500;
        allfi5 = cell2mat(fi5summary');
        hc = histcounts(allfi5(:,5),altbins) ./ size(fi5summary,2);
        plot(altbins(1:end-1),hc,'DisplayName',[parentfolder filesep fs(n).name]);
        title(sprintf('Collision per altitude (avg of %i runs)',size(fi5summary,2) ));
        xlabel('Altitude of collision (km)')
        ylabel('Count'); grid on;

    
    % Collision per altitude  (per year - imagesc)
    figure(1010 + n); clf;
        % altbins = 150:50:2500;
        allfi5 = cell2mat(fi5summary');
        hc2 = histcounts2(allfi5(:,5),allfi5(:,1),altbins,1:365/5:max(allfi5(:,1)));
        imagesc( (1:max(allfi5(:,1))) *5/365,altbins(1:end-1),hc2); 
        title('Collisions',fs(n).name,'interpreter','none'); colormap('jet'); 
            % frag5tmpvec = [n, iscata, p1_mass, p2_mass, (mean([norm(p1_r), norm(p2_r)])-radiusearthkm), size(debris1,1) + size(debris2,1), p1_radius, p2_radius, p1_objectclass, p2_objectclass, p1_cont, p2_cont];
        a = gca; a.YDir = 'normal'; xlabel('Years'); ylabel('Altitude (km)');
        colorbar;

    % Initial population per altitude
    figure(104); hold on;
        p1 = plot(paramSSEM.R02(1:end-1),tot_MC_mean(1,:),'DisplayName',[parentfolder filesep fs(n).name]);
        p2 = plot(paramSSEM.R02(1:end-1),tot_MC_mean(1,:) - mean(S_MC(1,:,:),3),'DisplayName','uncontrolled');
        p2.Color = p1.Color; p2.LineStyle = ':';
        title('Initial population');
        xlabel('Altitude (km)')
        ylabel('Count'); grid on;

    % Final population per altitude
    figure(105); hold on;
        % valid_tot_mean = tot_MC_mean;  % done earlier
        % valid_tot_mean(isnan(valid_tot_mean(:,1)),:) = []; % remove nan'd times
        valid_S = S_MC(nanIdx,:,:);
        % valid_S(isnan(valid_S(:,1)),:,:) = [];
        p1 = plot(paramSSEM.R02(1:end-1),valid_tot_mean(end,:),'DisplayName',[parentfolder filesep fs(n).name]);
        p2 = plot(paramSSEM.R02(1:end-1),valid_tot_mean(end,:) - mean(S_MC(end,:,:),3),'DisplayName','uncontrolled');
        p2.Color = p1.Color; p2.LineStyle = ':';
        title(sprintf('Final population (%i yrs)', floor(max(ts(~nanIdx)))));
        xlabel('Altitude (km)')
        ylabel('Count'); grid on;
    end
end

% legend('interpreter','none','fontsize',6)

end