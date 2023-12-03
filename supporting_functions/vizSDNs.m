function vizSDNs(S,D,Ns,B,paramSSEM,cfg)
% based on vizSDNB.m
% but modified for binned Ns  (size: nT, nShell, nBins, nMC)

% ts = [1:size(S,1)] * 5 / 365;  % yrs (main_mc)
if isfield(cfg,'saveMSnTimesteps')
    ts = [1:size(S,1)] * cfg.saveMSnTimesteps * 5 / 365;  % yrs
else
    ts = [1:size(S,1)] * 5 / 365;  % yrs
end

nshell = paramSSEM.N_shell;
nshellLabelInd = 1:10:nshell;
nshellLabel = cellstr(num2str(paramSSEM.R02(nshellLabelInd)'));

if isfield(paramSSEM,'NradiusEdges')
    % NmassEdges  = paramSSEM.NmassEdges;
    NmassEdges  = paramSSEM.NradiusEdges;  
else
    NmassEdges  = [0,inf];  
end

figure;  % for S and D (and B)
if ~isempty(B)
    subind = 3;
else
    subind = 2;
end

% per shell population via imagesc
subplot(2,subind,1);
imagesc(mean(S,3)'); 
a = gca; a.YDir = 'normal';
a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
title('S',sprintf('Avg of %i runs',size(S,3))); ylabel('Altitude (km)'); xlabel('Timesteps')

subplot(2,subind,2);
imagesc(mean(D,3)'); 
a = gca; a.YDir = 'normal';
a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
title('D',sprintf('Avg of %i runs',size(D,3))); ylabel('Altitude (km)'); xlabel('Timesteps')

% total number
subplot(2,subind,subind+1)
plot(ts,squeeze(sum(S,2))'); 
title('S'); ylabel('Number across shells'); xlabel('Years')

subplot(2,subind,subind+2);
plot(ts,squeeze(sum(D,2))');
title('D'); ylabel('Number across shells'); xlabel('Years')

if subind == 3  % plot RB data
    subplot(2,subind,subind);
    imagesc(mean(B,3)');
    a = gca; a.YDir = 'normal';
    a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
    title('RB',sprintf('Avg of %i runs',size(B,3))); ylabel('Altitude (km)'); xlabel('Timesteps')
    subplot(2,subind,subind*2);
    plot(ts,squeeze(sum(B,2))');
    title('RB'); ylabel('Number across shells'); xlabel('Years')
end

figure;  % for Ns
numNs = size(Ns,3);  % Ns:  nT x nShell x bins x nMC
for nind = 1:numNs
    N = squeeze(Ns(:,:,nind,:));
    subplot(2,numNs,nind);
    imagesc(mean(N,3)');
    a = gca; a.YDir = 'normal'; a.YTick = nshellLabelInd; a.YTickLabel = nshellLabel ;
    title('N',sprintf('%0.1e to %0.1e m',NmassEdges(nind),NmassEdges(nind+1))); 
    ylabel('Altitude (km)'); xlabel('Timesteps')

    subplot(2,numNs,nind+numNs);
    plot(ts,squeeze(sum(N,2))','Color',[0.2 0.5 0.9 0.2]);
    hold on;
    plot(ts,mean(squeeze(sum(N,2))',1),'LineWidth',2);
    title('N',sprintf('%0.1e to %0.1e m',NmassEdges(nind),NmassEdges(nind+1))); 
    ylabel('Number'); xlabel('Years')
    xlim([0,ts(end)])
end

figure;  % TOTAL
totnum = squeeze(sum(S,2)) + squeeze(sum(D,2)) + ...
     squeeze(sum(Ns,[2,3]));
if ~isempty(B)
    totnum = totnum + squeeze(sum(B,2));
end
plot(ts,totnum','Color',[0.2 0.5 0.9 0.2]);
hold on;
plot(ts,mean(totnum,2),'LineWidth',2);
title('Total number of objects'); 
ylabel('Number'); xlabel('Years')
xlim([0,ts(end)])

end