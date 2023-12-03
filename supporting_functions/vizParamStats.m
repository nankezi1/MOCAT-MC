function vizParamStats(param_mean, param_median, param_var, param_prc, paramSSEM)
% based on vizSDNB.m
% but modified for binned Ns  (size: nT, nShell, nBins, nMC)

% param_X: [nShell*nSpecies x 3:[bstar,mass,radius] x nT x nMC]

ts = [1:size(param_mean,3)] * 50 * 5 / 365;  % yrs
tsLabelInd = find(mod(ts,10) < 0.6);
tsLabel = cellstr(num2str(floor(ts(tsLabelInd)')));
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

%% EVERYTHING IN ONE LOOP?
for loop = 1:3
    switch loop
        case 1
            loopStr = 'B^*';
        case 2
            loopStr = 'm';
        case 3
            loopStr = 'r';
    end
    for ind = 1:nSpecies
        if ind == 1 | ind == 3
            figure;
        end
        if ind < 3  % for S/D
            continue; % skip for now
            subplot(3,2,ind); % mean
            imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca;
            if loop == 1
                a.CLim = prctile(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];  % B* colormap is 5th-95th percentile
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            % title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)));
            ttext = ['$\bar ' loopStr '_' speciesNames{ind} '$'];
            title(ttext, '', 'interpreter','latex');
            ylabel('Altitude (km)'); %xlabel('Years')
            colorbar;
            subplot(3,2,2+ind); % median
            imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca; 
            if loop == 1
                a.CLim = prctile(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4)));
            ttext = ['$median(' loopStr '_' speciesNames{ind} ')$'];
            title(ttext, '', 'interpreter','latex');
            ylabel('Altitude (km)'); %xlabel('Years')
            colorbar;
            subplot(3,2,4+ind); % var
            imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca;
            if loop == 1
                a.CLim = prctile(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)));
            ttext = ['$Var(' loopStr '_' speciesNames{ind} ')$'];
            title(ttext, '', 'interpreter','latex');
            ylabel('Altitude (km)'); xlabel('Years')
            colorbar;
        else  % Ns
            subplot(3,nSpecies-2,ind-2); % mean
            imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca; 
            if loop == 1
                a.CLim = prctile(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];  % colormap is 5th-95th percentile
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)),'interpreter','none');
            ttextSplit = strsplit( speciesNames{ind});
            ttextRng = strsplit(ttextSplit{2},'_');
            ttext = ['$\bar ' loopStr '_N$'];
            title(ttext, ['N: ' ttextRng{1} ' to ' ttextRng{2}], 'interpreter','latex');
            ylabel('Altitude (km)'); xlabel('Years')
            colorbar;
            
            subplot(3,nSpecies-2,nSpecies+ind-4); % median
            imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca; ...
            if loop == 1 ...
                a.CLim = prctile(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4)),'interpreter','none');
            ttext = ['$median(' loopStr '_N)$'];
            title(ttext,  'interpreter','latex');
            ylabel('Altitude (km)'); xlabel('Years')
            colorbar;
            
            subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % var
            imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
            a = gca; 
            if loop == 1
                a.CLim = prctile(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6];
            end
            a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
            a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
            %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none');
            ttext = ['$Var(' loopStr '_N)$'];
            title(ttext, 'interpreter','latex');
            ylabel('Altitude (km)'); xlabel('Years')
            colorbar;

            % percentiles (prc) defined: MC2SSEM_population_dist_binned.m
            if ~isempty(param_prc)
                q05 = cellfun(@(x) double(x(1)), param_prc);
                q25 = cellfun(@(x) double(x(2)), param_prc);
                q75 = cellfun(@(x) double(x(3)), param_prc);
                q95 = cellfun(@(x) double(x(4)), param_prc);
                iqr2575 = q75-q25;
                iqr0595 = q95-q05;
    
                subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % IQR 25-75
                imagesc(squeeze(mean(iqr2575(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
                a = gca; 
                a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
                a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
                %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none');
                ttext = ['$IQR-25-75(' loopStr '_N)$'];
                title(ttext, 'interpreter','latex');
                ylabel('Altitude (km)'); xlabel('Years')
                colorbar;
    
                % IQR: 5-95
                subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % IQR 25-75
                imagesc(squeeze(mean(iqr0595(1+(ind-1)*nshell:nshell*ind,loop,:,:),4,'omitnan')));
                a = gca; 
                a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
                a.XTick = tsLabelInd; a.XTickLabel = tsLabel;
                %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none');
                ttext = ['$IQR-05-95(' loopStr '_N)$'];
                title(ttext, 'interpreter','latex');
                ylabel('Altitude (km)'); xlabel('Years')
                colorbar;
            end
        end
    end
end

%%
return;


% % BSTAR
% for ind = 1:nSpecies
%     if ind == 1 | ind == 3
%         figure;
%     end
%     if ind < 3  % B* for S
%         subplot(3,2,ind); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);  % colormap is 5th-95th percentile
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         % title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)));
%         ttext = ['$\bar B^*_' speciesNames{ind} '$'];
%         title(ttext, '', 'interpreter','latex'); 
%         ylabel('Altitude (km)'); %xlabel('Years')
%         colorbar;
%         subplot(3,2,2+ind); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4))); 
%         ttext = ['$median(B^*_' speciesNames{ind} ')$'];
%         title(ttext, '', 'interpreter','latex'); 
%         ylabel('Altitude (km)'); %xlabel('Years')
%         colorbar;
%         subplot(3,2,4+ind); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);        
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4))); 
%         ttext = ['$Var(B^*_' speciesNames{ind} ')$'];
%         title(ttext, '', 'interpreter','latex'); 
%         ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     else  % B* for D
%         subplot(3,nSpecies-2,ind-2); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_mean(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);  % colormap is 5th-95th percentile
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         %title([speciesNames{ind} ' mean B*'],sprintf('Avg of %i runs',size(param_mean,4)),'interpreter','none'); 
%         ttextSplit = strsplit( speciesNames{ind});
%         ttextRng = strsplit(ttextSplit{2},'_');
%         ttext = ['$\bar B^*_N$'];
%         title(ttext, ['N: ' ttextRng{1} ' to ' ttextRng{2}], 'interpreter','latex'); 
%         ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,nSpecies+ind-4); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_median(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         %title([speciesNames{ind} ' median B*'],sprintf('Avg of %i runs',size(param_median,4)),'interpreter','none'); 
%         ttext = ['$median(B^*_N)$'];
%         title(ttext,  'interpreter','latex'); 
%         ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),4,'omitnan')),...
%             prctile(param_var(1+(ind-1)*nshell:nshell*ind,1,:,:),[5,95],'all')'+[0,1e-6]);
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         %title([speciesNames{ind} ' var B*'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none'); 
%         ttext = ['$Var(B^*_N)$'];
%         title(ttext, 'interpreter','latex'); 
%         ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     end
% end
% 
% 
% % MASS
% for ind = 1:nSpecies
%     if ind == 1 | ind == 3
%         figure;
%     end
%     if ind < 3  % B* for S
%         subplot(3,2,ind); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' mean mass'],sprintf('Avg of %i runs',size(param_mean,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,2,2+ind); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' median mass'],sprintf('Avg of %i runs',size(param_median,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,2,4+ind); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' var mass'],sprintf('Avg of %i runs',size(param_var,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     else  % B* for D
%         subplot(3,nSpecies-2,ind-2); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' mean mass'],sprintf('Avg of %i runs',size(param_mean,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,nSpecies+ind-4); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' median mass'],sprintf('Avg of %i runs',size(param_median,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,2,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' var mass'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     end
% end
% 
% 
% % RADIUS
% for ind = 1:nSpecies
%     if ind == 1 | ind == 3
%         figure;
%     end
%     if ind < 3  % B* for S
%         subplot(3,2,ind); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         title([speciesNames{ind} ' mean radius'],sprintf('Avg of %i runs',size(param_mean,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,2,2+ind); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' median radius'],sprintf('Avg of %i runs',size(param_median,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,2,4+ind); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' var radius'],sprintf('Avg of %i runs',size(param_var,4))); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     else  % B* for D
%         subplot(3,nSpecies-2,ind-2); % mean
%         imagesc(squeeze(mean(param_mean(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' mean radius'],sprintf('Avg of %i runs',size(param_mean,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,nSpecies+ind-4); % median
%         imagesc(squeeze(mean(param_median(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' median radius'],sprintf('Avg of %i runs',size(param_median,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%         subplot(3,nSpecies-2,2*(nSpecies-2)+ind-2); % var
%         imagesc(squeeze(mean(param_var(1+(ind-1)*nshell:nshell*ind,3,:,:),4,'omitnan')));
%         a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
%         a.XTick = tsLabelInd; a.XTickLabel = tsLabel; 
%         title([speciesNames{ind} ' var radius'],sprintf('Avg of %i runs',size(param_var,4)),'interpreter','none'); ylabel('Altitude (km)'); xlabel('Years')
%         colorbar;
%     end
% end

end