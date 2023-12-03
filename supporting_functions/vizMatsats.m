function vizMatsats(matsat, varargin)
% Visualize given matsat into shell vs SSEM types

% if isfield(cfg,'paramSSEM') && isfield(cfg.paramSSEM,'R02')
%     R02 = cfg.paramSSEM.R02;
%     re = cfg.paramSSEM.re;
% else
    R02 = 200:50:2000;
    re = 6378.137;
% end

if nargin > 1
    figure(varargin{1})
else
    figure;
end

idx_a = 1; idx_ecco = 2; idx_inclo = 3;
idx_nodeo = 4; idx_argpo = 5; idx_mo = 6;
idx_bstar = 7; idx_mass = 8; idx_radius = 9;
idx_error = 10; idx_controlled = 11; idx_a_desired = 12;   % [0,0,0]
idx_missionlife = 13; idx_constel = 14; idx_date_created = 15;  % [0,0,nan]
idx_launch_date = 16; idx_r = [17 18 19]; idx_v = [20 21 22]; % date: JD,  RV:0
idx_objectclass = 23; idx_ID = 24;


[g1,g2,g3] = getZeroGroups(matsat);
idxCtrl = find(matsat(:,idx_controlled) == 1);
idxNoCtrl = find(matsat(:,idx_controlled) == 0);
Scount = histcounts(matsat(intersect(g1.allclass, idxCtrl), idx_a) , 1 + R02 / re);
Dcount = histcounts(matsat(intersect(g1.allclass, idxNoCtrl), idx_a) ,  1 + R02 / re);
Bcount = histcounts(matsat(g2.allclass, idx_a) ,  1 + R02 / re);
Ncount = histcounts(matsat(g3.allclass, idx_a) ,  1 + R02 / re);

%% FIGURE: Initial population from matsat = cfg.mat_sats
% figure;
% plot(R02(1:end-1) + diff(R02), Scount, '--x'); hold on;
% plot(R02(1:end-1) + diff(R02), Dcount, '--x');
% plot(R02(1:end-1) + diff(R02), Ncount, '--x');
% plot(R02(1:end-1) + diff(R02), Bcount, '--x');
% plot(R02(1:end-1) + diff(R02), Scount + Dcount + Ncount + Bcount, 'k--o');
% legend('S','D','N','B','Total')
% 
% xlabel('Shell Altitude (km)')
% ylabel('Count')
% % title('Initial Population',sprintf('Epoch: %s', cfg.time0))
% ylim([0,max(Scount + Dcount + Ncount + Bcount) * 1.2]);

fprintf('Mean count per shell: S %i, D %i, N %i, RB %i \n', round(mean([Scount; Dcount; Ncount; Bcount],2)))


%% FIGURE: SUBPLOT HISTOGRAM
c = lines(3); 

% figure;
subplot(311);
bar(R02(1:end-1)+ diff(R02), Scount,'facealpha',0.7,'FaceColor',c(1,:)); title('S')
subplot(312);
bar(R02(1:end-1)+ diff(R02), Dcount,'facealpha',0.7,'FaceColor',c(2,:)); title('D')
subplot(313);
bar(R02(1:end-1)+ diff(R02), Ncount,'facealpha',0.7,'FaceColor',c(3,:)); title('N')
xlabel('Alt (km)')

% figure;  % STACKED
% subplot(211);
% bar(paramSSEM.R02(1:end-1), SS(end,:),'facealpha',0.4,'FaceColor',c(1,:)); hold on;
% bar(paramSSEM.R02(1:end-1), DD(end,:),'facealpha',0.5,'FaceColor',c(2,:)); 
% bar(paramSSEM.R02(1:end-1), NN(end,:),'facealpha',0.6,'FaceColor',c(3,:)); 
% legend('S','D','N')
% title('S,D,N per shell at t=tf (overlapping)')
% subplot(212);
% bar(paramSSEM.R02(1:end-1), [SS(end,:); DD(end,:); NN(end,:)]','stacked');
% legend('S','D','N');
% title('S,D,N per shell at t=tf (stacked)')

end



function [g1,g2,g3] = getZeroGroups(inmatsat)
idx_mass = 8; idx_radius = 9; idx_objectclass = 23;

% mass vs radius
g1.allclass = []; % group 1: payloads; logical index
g2.allclass = []; % group 2: RBs
g3.allclass = []; % group 3: all debris

for ii = 1:12
    msinds = inmatsat(:,idx_objectclass) == ii; % all obj w/ objclass
    if ii == 1
        g1.allclass = find(msinds);                                 % all payload entries
        g1.zr = g1.allclass(inmatsat(g1.allclass,idx_radius) == 0); % index of zero radius
        g1.zm = g1.allclass(inmatsat(g1.allclass,idx_mass) == 0);   % index of zero mass
        g1.nz = setdiff(g1.allclass, union(g1.zr,g1.zm));           % index of non-zero radius and mass
        g1.nzno = g1.nz(~isoutlier(inmatsat(g1.nz,idx_radius)) ...  % non-zero, non-outlier
            & ~isoutlier(inmatsat(g1.nz,idx_mass)) );
    elseif ii == 5 | ii == 6
        g2.allclass = find(msinds);
        g2.zr = g2.allclass(inmatsat(g2.allclass,idx_radius) == 0); % all RB entries
        g2.zm = g2.allclass(inmatsat(g2.allclass,idx_mass) == 0);
        g2.nz = setdiff(g2.allclass, union(g2.zr,g2.zm));
        g2.nzno = g2.nz(~isoutlier(inmatsat(g2.nz,idx_radius)) ...
            & ~isoutlier(inmatsat(g2.nz,idx_mass)) );
    else
        g3.allclass = union(g3.allclass, find(msinds));             % all debris entries
    end
end
g3.zr = g3.allclass(inmatsat(g3.allclass,idx_radius) == 0);
g3.zm = g3.allclass(inmatsat(g3.allclass,idx_mass) == 0);
g3.nz = setdiff(g3.allclass, union(g3.zr, g3.zm));
g3.nzno = g3.nz(~isoutlier(inmatsat(g3.nz,idx_radius)) ...
    & ~isoutlier(inmatsat(g3.nz,idx_mass)) );
end