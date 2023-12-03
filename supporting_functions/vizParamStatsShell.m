function vizParamStatsShell(param_mean, param_median, param_var, paramSSEM, shellnums)
% based on vizSDNB.m, vizParamStats
% for specific shells in shellnums (2D graph instead of imagesc)
% also, allow for empty param_X, to skip plotting it

% param_X: [nShell*nSpecies x 3:[bstar,mass,radius] x nT x nMC]

paramsOn = [~isempty(param_mean), ~isempty(param_median), ~isempty(param_var)];

if sum(paramsOn) < 1
    warning('No param_X given')
    return;
end

ts = [1:size(param_mean,3)] * 60 * 5 / 365;  % yrs   % MCconfig.saveMSnTimesteps, dt_days
% tsLabelInd = find(mod(ts,10) < 0.6);
% tsLabel = cellstr(num2str(floor(ts(tsLabelInd)')));
nshell = paramSSEM.N_shell;
nshellLabelInd = 1:10:nshell;
nshellLabel = cellstr(num2str(paramSSEM.R02(nshellLabelInd)'));
nSpecies = size(param_mean,1) / nshell ;

if ~isempty(paramSSEM.NradiusEdges)
    Nbinedges = paramSSEM.NradiusEdges;
else
    Nbinedges = paramSSEM.NmassEdges;
end

speciesNames = {'S','D'};
for ind = 3:nSpecies
    speciesNames{ind} = sprintf('N %0.1e_%0.1e', Nbinedges(ind-2:ind-1));
end

%% EVERYTHING IN ONE LOOP!

for loop = 1:3
    switch loop
        case 1
            loopStr = 'B^*';
            loopUnit = '(1/R_E)';
        case 2
            loopStr = 'm';
            loopUnit = '(kg)';
        case 3
            loopStr = 'r';
            loopUnit = '(m)';
    end
    for ind = 1:nSpecies
        if ind == 3  % |  ind == 1 
            figure;
        end
        if ind < 3  % B* for S
            continue; % skip for now

            if paramsOn(1)
            subplot(3,2,ind); % mean
            plot(ts,(squeeze(param_mean((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (means) individually
            % plot(ts,(squeeze(mean(param_mean((ind-1)*nshell+shellnums,loop,:,:),4,'omitnan')))); % mean of runs (means)
            a = gca;
%             a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            % title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)));
            ttext = ['$\bar ' loopStr '_' speciesNames{ind} '$'];
            title(ttext, '', 'interpreter','latex');
            ylabel([ttext(1:end-1), '\ \ ' loopUnit '$'], 'interpreter','latex'); %xlabel('Years')
            end
            
            if paramsOn(2)
            subplot(3,2,2+ind); % median
            plot(ts,(squeeze(param_median((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (median) individually
            a = gca; 
%             a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4)));
            ttext = ['$median(' loopStr '_' speciesNames{ind} ')$'];
            title(ttext, '', 'interpreter','latex');
            ylabel([ttext(1:end-1), '\ \ ' loopUnit '$'], 'interpreter','latex'); %xlabel('Years')
            end
            if paramsOn(3)
            subplot(3,2,4+ind); % var
            plot(ts,(squeeze(param_var((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (var) individually
            a = gca;
            %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)));
            ttext = ['$Var(' loopStr '_' speciesNames{ind} ')$'];
            title(ttext, '', 'interpreter','latex');
            ylabel(ttext, 'interpreter','latex'); xlabel('Years')
            end
        else  % B* for D
            if paramsOn(1)
            subplot(3,nSpecies-2,ind-2); % mean
            plot(ts,(squeeze(param_mean((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (mean) individually
            a = gca; 
%             a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)),'interpreter','none');
            ttextSplit = strsplit( speciesNames{ind});
            ttextRng = strsplit(ttextSplit{2},'_');
            ttext = ['$\bar ' loopStr '_N$'];
            title(ttext, ['N: ' ttextRng{1} ' to ' ttextRng{2}], 'interpreter','latex');
            ylabel([ttext(1:end-1), '\ \ ' loopUnit '$'], 'interpreter','latex'); %xlabel('Years')
            end
            if paramsOn(2)
            subplot(3,nSpecies-2,nSpecies+ind-4); % median
            plot(ts,(squeeze(param_median((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (median) individually
            a = gca;
            %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4)),'interpreter','none');
            ttext = ['$median(' loopStr '_N)$'];
            title(ttext,  'interpreter','latex');
            ylabel([ttext(1:end-1), '\ \ ' loopUnit '$'], 'interpreter','latex'); %xlabel('Years')
            end
            if paramsOn(3)
            subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % var
            plot(ts,(squeeze(param_var((ind-1)*nshell+shellnums,loop,:,:)))); % all runs (var) individually
            a = gca; 

%             a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%             a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none');
            ttext = ['$Var(' loopStr '_N)$'];
            title(ttext, 'interpreter','latex');
            ylabel(ttext, 'interpreter','latex'); xlabel('Years')
            end
        end
    end
end

end