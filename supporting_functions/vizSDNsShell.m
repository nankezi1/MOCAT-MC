function vizSDNsShell(S,D,Ns,paramSSEM,shellnums)
% based on vizSDNB.m, vizSDNs
% but modified for binned Ns  (size: nT, nShell, nBins, nMC)
% and for specific shell (2D graph)

% ts = [1:size(S,1)] * 5 / 365;  % yrs (main_mc)
ts = [1:size(S,1)] * 60 * 5 / 365;  % yrs (main_mc2; see MCconfig.saveMSnTimesteps and .dt_days)

% nshell = paramSSEM.N_shell;
% nshellLabelInd = 1:10:nshell;
% nshellLabel = cellstr(num2str(paramSSEM.R02(nshellLabelInd)'));

% NmassEdges  = paramSSEM.NmassEdges;
NmassEdges  = paramSSEM.NradiusEdges;  

figure;  % for S and D

S = S(:,shellnums,:);
D = D(:,shellnums,:);
Ns = Ns(:,shellnums,:,:);
j = jet(shellnums);

% per shell population via imagesc
subplot(221);
% imagesc(mean(S,3)'); 
plot(ts,mean(S,3)'); colororder(hsv(size(S,2)));
title('S',sprintf('Avg of %i runs, per shell',size(S,3))); 
legend(strcat(num2str(paramSSEM.R02(shellnums)'), ' km'));
ylabel('Population'); 
xlabel('Timesteps')

subplot(222); cla
% imagesc(mean(D,3)'); 
plot(ts,mean(D,3)'); colororder(jet(size(D,2)));

    for ind = 1:size(D,2)
        hold on;
        plot(ts,squeeze(D(:,ind,:)),':','Color',j(ind,:));
    end
title('D',sprintf('Avg of %i runs, per shell',size(D,3))); 
% ylabel('Altitude (km)'); 
xlabel('Timesteps')

% total number
subplot(223)
plot(ts,squeeze(sum(S,2))'); 
title('S'); ylabel('Total population per MC'); 
xlabel('Years')

subplot(224);
plot(ts,squeeze(sum(D,2))');
title('D'); ylabel('Total population per MC');
xlabel('Years')


figure;  % for Ns
numNs = size(Ns,3);  % Ns:  nT x nShell x bins x nMC
for nind = 1:numNs
    N = squeeze(Ns(:,:,nind,:));
    subplot(2,numNs,nind);
%     imagesc(mean(N,3)');
    plot(ts,mean(N,3)'); colororder(jet(size(N,2)));
    title('N',sprintf('%0.1e to %0.1e m',NmassEdges(nind),NmassEdges(nind+1))); 
    ylabel('Population'); xlabel('Timesteps')
        
        for ind = 1:size(D,2)
            hold on;
            plot(ts,squeeze(D(:,ind,:)),':','Color',j(ind,:));
        end

    subplot(2,numNs,nind+numNs);
    plot(ts,squeeze(sum(N,2))');
    title('N',sprintf('%0.1e to %0.1e m',NmassEdges(nind),NmassEdges(nind+1))); 
    ylabel('Total population'); xlabel('Years')
end

end