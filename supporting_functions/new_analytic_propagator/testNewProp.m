close all; 
clear all;
clc;

global tumin radiusearthkm xke j2 j3 j4 j3oj2


% Test script to check SGP4 (spg4_ecf.m) behavior
% addpath('tle_plots')
addpath(genpath(pwd))

% if isfile('sats.mat')
%     load sats.mat
% else
    name='all_2022.txt';
    [sats,~,~]=load_all_tles_3lines(name);
% end

% [sats]=add_mass_radius(sats);

% compute all Perigees to filter non-LEO objects
n_sats=length(sats);
for k=1:n_sats
    a_all(k,1)=sats{k}.a*radiusearthkm;
    ap_all(k,1)=(sats{k}.a*radiusearthkm)*(1-sats{k}.ecco);
    aa_all(k,1)=(sats{k}.a*radiusearthkm)*(1+sats{k}.ecco);
    % above is equivalent to:
    % aa_all(k,1)=(sats{k}.alta*radiusearthkm);
    % ap_all(k,1)=(sats{k}.altp*radiusearthkm);
end


% Choose low perigee satellite; 
lowsat = sats((ap_all > 300+radiusearthkm) & (ap_all < 350+radiusearthkm) ...
    & (aa_all > 500+radiusearthkm) & (aa_all < 1000+radiusearthkm));  

lowsat = lowsat{1};   % ind 5938: ap: 328km; per: 986km

time0=datetime(lowsat.jdsatepoch,'convertfrom','juliandate');

dt = 24*60;                            % hrs
tsince = 0 : dt : 24*60*100;    % for 10 years
n_time=length(tsince);

X_ecf = zeros(6,n_time); % ECEF via sgp4_ecf.mRi
X_eci = zeros(6,n_time); % ECI via sgp4.m

sat1=lowsat;
[sat1]=int_sgp4_mit(sat1);
for n=1:n_time
   [tmpsat, X_ecf_temp, ~,~] = spg4_ecf(lowsat,tsince(n));           % ECF
   [~, X_eci_temp_r, X_eci_temp_v] = sgp4(lowsat,tsince(n));    % ECI
  %[sat1,X_eci_temp_r2,X_eci_temp_v2,oe_mean_out(:,n)]=sgp4_mit(sat1,dt);
  [sat1,X_eci_temp_r2,X_eci_temp_v2]=spg4_mitv2(sat1,tsince(n));

   if size(X_ecf_temp,1) ~= 6
       warning('SGP4_ecf error; output size is not 6x1');
       fprintf('satrec error code: %i \n',tmpsat.error);
       break;
   else
       X_ecf(:,n) = X_ecf_temp;
       X_eci(:,n) = [X_eci_temp_r'; X_eci_temp_v'];
       X_eci2(:,n) = [X_eci_temp_r2; X_eci_temp_v2];

       if mod(n,round(n_time/10)) == 0
           fprintf('%0.1f%%,  satrec error: %i \n', n/n_time*100, tmpsat.error)
       end
   end
end

r_ecf = sqrt(X_ecf(1,:).^2 + X_ecf(2,:).^2 + X_ecf(3,:).^2);
r_eci = sqrt(X_eci(1,:).^2 + X_eci(2,:).^2 + X_eci(3,:).^2);
r_eci2 = sqrt(X_eci2(1,:).^2 + X_eci2(2,:).^2 + X_eci2(3,:).^2);


figure(1); clf;
subplot(3,1,1)
plot(tsince/60/24, r_ecf - radiusearthkm, 'r-'); 
legend('r$_{ecf}$'); ylabel('r (km)')
subplot(3,1,2);
plot(tsince/60/24, r_eci - radiusearthkm, 'b-'); 
legend('r$_{eci}$'); ylabel('r (km)')
subplot(3,1,3);
plot(tsince/60/24, r_eci2 - radiusearthkm, 'k-');
legend('r$_{eci,mit}$'); xlabel('days'); ylabel('r (km)')
%ylim([0,lowsat.alta * radiusearthkm * 1.2])


figure(10); clf;
plot(tsince/60/24, r_ecf - radiusearthkm, 'r.-'); hold on;
plot(tsince/60/24, r_eci - radiusearthkm, 'b.-');
plot(tsince/60/24, r_eci2 - radiusearthkm, 'k.-');
xlabel('days'); ylabel('r (km)')
legend('r$_{ecf}$','r$_{eci}$','r$_{eci,mit}$')
%ylim([0,lowsat.alta * radiusearthkm * 1.2])
