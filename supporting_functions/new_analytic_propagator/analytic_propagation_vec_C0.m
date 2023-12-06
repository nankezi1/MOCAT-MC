function [out_oe,errors]=analytic_propagation_vec(input_oe,param)

errors=zeros(size(input_oe,1),1);
% this code includes the solution for the mean elements
% as a function of time from Ref. 1.
%
% [Ref. 1] Martinusi, Vladimir, Lamberto Dell Elce, and Gaëtan Kerschen.
% "Analytic propagation of near-circular satellite orbits
% in the atmosphere of an oblate planet." Celestial Mechanics
% and Dynamical Astronomy 123, no. 1 (2015): 85-103.

re = param.req;
J2 = param.j2;
mu = param.mu;

% if statement ensure we get reasonable values for rho. Rho=0 break the
% model below. Also altitude greater than 2300 km breaks the model.
a_0=input_oe(:,1);
a_minus_re = a_0-re;
rho_0 = zeros(length(a_0),1);

% if strcmpi(param.density_profile,'JB2008')
%     
%     rho_0 = zeros(length(a_0),1);
%     check_above = a_minus_re>param.alt(end,1);
%     check_below = a_minus_re<param.alt(1,1);
%     check_in_range = ~check_above & ~check_below;
% 
%     % if statement ensure we get reasonable values for rho. Rho=0 break the
%     % model below. Also altitude greater than 2300 km breaks the model.
% %     if (input_oe(1,1)-re)<2000 && (input_oe(1,1)-re)>200% changes from 200-1100
%     
%         %rho_0 = interp2(param.dens_times,param.alt,param.dens_value,param.jd,input_oe(1,1)-re)*(1000)^3;
%     %     rho_0 =qinterp2((param.alt)',(param.dens_times)',param.dens_value',input_oe(1,1)-re, param.jd,1)*(1000)^3;
%         %     F107=param.F107;
%         %     Ap=param.Ap;
%         %     % density model
%         %     h=input_oe(1,1);
%         %     T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
%         %     m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
%         %     H = T / m;
%         %     rho_0 = 6e-10* exp ( - ( input_oe(1,1)-re- 175 ) / H ) *(1000)^3;%units kg/km^3
%     %     rho_old = qinterp2((param.alt)',(param.dens_times)',param.dens_value',input_oe(1,1)-re, param.jd,1)*(1000)^3;
%     %     rho_grid =param.dens_grid(input_oe(1,1)-re, param.jd);
%     rho_0(check_in_range) = lininterp2_vec_v2(param.alt(:,1),param.dens_times(1,:),param.dens_value,a_minus_re(check_in_range), param.jd)*1e9;
%         % https://www.mathworks.com/matlabcentral/fileexchange/28376-faster-linear-interpolation
%     
%     %     fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
% %     else
% %         if (input_oe(1,1)-re)>1100
%             %         aa=1100+re;
%             %         F107=param.F107;
%             %         Ap=param.Ap;
%             %         % density model
%             %         h=aa;
%             %         T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
%             %         m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
%             %         H = T / m;
%             %         rho_0 = 6e-10* exp ( - ( aa-re- 175 ) / H ) *(1000)^3;%units kg/km^3
%             % rho_0 = interp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
%     rho_0(check_above) = lininterp1_vec(param.dens_times(1,:),param.dens_value(end,:),param.jd)*1e9;
%     %         rho_old = interp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
%     %         fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
%     
% %         else
%             %         aa=200+re;
%             %         F107=param.F107;
%             %         Ap=param.Ap;
%             %         % density model
%             %         h=aa;
%             %         T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
%             %         m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
%             %         H = T / m;
%             %         rho_0 = 6e-10* exp ( - ( aa-re- 175 ) / H ) *(1000)^3;%units kg/km^3
%             % rho_0 = interp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
%     rho_0(check_below) = lininterp1_vec(param.dens_times(1,:),param.dens_value(1,:),param.jd)*1e9;
%     %         rho_old = interp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
%     %         fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
% %         end
% %     end
% 
% elseif strcmpi(param.density_profile,'static')
%     rho_0 = densityexp_vec(a_minus_re)*1e9;
% end

%constants
%C_0=param.C_0*rho_0;%1/2*Cd*s_ref/m*rho_0

%check these unit unit of 1/re and using 1/(6378.137 km) to coтvert to km
%units

% rho_reference for Bstar is 0.157 in units of kg/(m^2 * re).
% Using (1000^2 m^2/1 km^2) to convert A from m^2 to km^2, and dividing by 0.157 to scale
% according to newly computed density rho_0 (kg/km^3)% 1/(6378.137 km)*(1000^2 m/1 km) to conver to kg/km^3

% C_0 = (param.Bstar/(1e6*0.157)).*rho_0;
% min_C = min(C_0)
C_0 = logspace(-8,-22,numel(rho_0))';

threek2 = 3*mu*J2*re^2/2;

t = param.t;
t_0 = param.t_0;

%initial conditions
e_0 = input_oe(:,2);
inc_0 = input_oe(:,3);
bigO_0 = input_oe(:,4);
omega_0 = input_oe(:,5);
Mo_0 = input_oe(:,6);

c = cos(inc_0);
c_sq = c.^2;

n_0 = sqrt(mu)*a_0.^(-3/2);

alpha0_sq = (e_0./sqrt(a_0)).^2;

beta_0 = (sqrt(3)/2)*e_0;

beta0_sq = beta_0.^2;

n0_a0_dt = n_0.*a_0*(t-t_0);

tan_atan_beta0 = max(tan(atan(beta_0)-beta_0.*n0_a0_dt.*C_0),0); %place lower limit on eccentricity and semi-major axis reduction
a = (a_0./beta0_sq).*tan_atan_beta0.^2;
e = (2/sqrt(3))*tan_atan_beta0;     

% p_0=a_0.*(1-e_0.^2);

% Cu=(3/4)*n_0*J2.*(re./p_0).^2;

%%%%%% No exception for small C_0
%Save some variables to avoid repetition of operations
a_sq = a.^2;
four_thirds_over_a_cb = 4/3./(a_sq.*a);
a0_sq = a_0.^2;
four_thirds_over_a0_cb = 4/3./(a0_sq.*a_0);
alpha0sq_over_asq = alpha0_sq./a_sq;
alpha0sq_over_a0sq = alpha0_sq./a0_sq;

Mo1 = 0.125*(4./a - 4./a_0 + 3*alpha0_sq.*log(a./a_0))./C_0 + threek2/(16*mu)*(3*c_sq-1).*(1.5*(alpha0sq_over_asq-alpha0sq_over_a0sq) + four_thirds_over_a_cb-four_thirds_over_a0_cb)./C_0 + Mo_0;

five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 = (2.5*(alpha0sq_over_asq-alpha0sq_over_a0sq) + four_thirds_over_a_cb-four_thirds_over_a0_cb)./C_0;

omega1 = threek2/(16*mu)*(5*c_sq-1).*five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 + omega_0;

bigO1 = -threek2/(8*mu)*c.*five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 + bigO_0;
%%%%%%%
figure(11)
plot(C_0,Mo1,'bo')

figure(12)
plot(C_0,omega1,'bo')

figure(13)
plot(C_0,bigO1,'bo')

%Save some variables to avoid repetition of operations
check_C0 = C_0<1e-18;
find_smallC0 = find(check_C0);
find_largeC0 = find(~check_C0);

%Equations for large C_0
a_large = a(find_largeC0);
a0_large = a_0(find_largeC0);

a_sq_large = a_large.^2;
a0_sq_large = a0_large.^2;
a_cb_large = a_sq_large.*a_large;
a0_cb_large = a0_sq_large.*a0_large;
alpha0_sq_large = alpha0_sq(find_largeC0);
four_thirds_over_a_cb_large = 4/3./a_cb_large;
four_thirds_over_a0_cb_large = 4/3./a0_cb_large;
alpha0sq_over_asq_large = alpha0_sq_large./a_sq_large;
alpha0sq_over_a0sq_large = alpha0_sq_large./a0_sq_large;

C0_large = C_0(find_largeC0);

Mo_first_term_large = (4./a_large - 4./a0_large + 3*alpha0_sq_large.*log(a_large./a0_large))./C0_large;
Mo_second_term_large = (1.5*(alpha0sq_over_asq_large - alpha0sq_over_a0sq_large) + four_thirds_over_a_cb_large - four_thirds_over_a0_cb_large)./C0_large;

omega_term_large = (2.5*(alpha0sq_over_asq_large - alpha0sq_over_a0sq_large) + four_thirds_over_a_cb_large - four_thirds_over_a0_cb_large)./C0_large;                                 

%Equations for small C_0
a_small = a(find_smallC0);
a0_small = a_0(find_smallC0);

a2_small = a_small.^2;
a02_small = a0_small.^2;
a3_small = a2_small.*a_small;
a03_small = a02_small.*a0_small;
a_a0 = a_small.*a0_small;
a2_a02 = a2_small.*a02_small;
a3_a03 = a3_small.*a03_small;

beta0_sq_small = beta0_sq(find_smallC0);
beta0_4 = beta0_sq_small.*beta0_sq_small;

d1 = -2*ones(numel(find_smallC0),1);
d2 = 1+3*beta0_sq_small;
d3 = -4/3*beta0_sq_small.*(2+3*beta0_sq_small);
d4 = 1/3*beta0_sq_small.*(15*beta0_4+15*beta0_sq_small+2);
d5 = -2/15*beta0_sq_small.*(45*beta0_4+60*beta0_sq_small+17);
d_vec = [d1,d2,d3,d4,d5];

g1 = d1;
g2 = -1+beta0_sq_small;
g3 = -2/3*(1+beta0_4);
g4 = 1/6*(1-beta0_sq_small).*(3*beta0_4+4*beta0_sq_small+3);
g5 = -2/15*beta0_sq_small.*(3*beta0_4.*beta0_4+2*beta0_4.*beta0_sq_small+2*beta0_sq_small+3);
g_vec = [g1,g2,g3,g4,g5];

C0_small = C_0(find_smallC0);
C02_small = C0_small.^2;
C03_small = C02_small.*C0_small;
C04_small = C02_small.*C02_small;

n0_a0_dtsmall = n0_a0_dt(find_smallC0);
n0_a0_dt2 = n0_a0_dtsmall.*n0_a0_dtsmall;
n0_a0_dt3 = n0_a0_dt2.*n0_a0_dtsmall;
n0_a0_dt4 = n0_a0_dt2.*n0_a0_dt2;
n0_a0_dt5 = n0_a0_dt4.*n0_a0_dtsmall;

n0_a0_dt_C0vec = [n0_a0_dtsmall,n0_a0_dt2.*C0_small,n0_a0_dt3.*C02_small,n0_a0_dt4.*C03_small,n0_a0_dt5.*C04_small];

one_plus_beta0_sq = 1+beta0_sq_small;

one_over_C0_a0_minus_a = -a0_small.*one_plus_beta0_sq.*sum(n0_a0_dt_C0vec.*d_vec,2);

one_over_C0_one_over_a_minus_one_over_a0 = one_over_C0_a0_minus_a./a_a0;
one_over_C0_one_over_a2_minus_one_over_a02 = one_over_C0_a0_minus_a./a2_a02.*(a0_small+a_small);
one_over_C0_one_over_a3_minus_one_over_a03 = one_over_C0_a0_minus_a./a3_a03.*(a02_small+a_a0+a2_small);

one_over_C0_log_a_over_a0 = n0_a0_dtsmall.*one_plus_beta0_sq.*sum(n0_a0_dt_C0vec.*g_vec,2).*C0_small;

alpha0_sq_small = alpha0_sq(find_smallC0);

one_over_C0_alpha0sq_over_asq_a0sq = alpha0_sq_small.*one_over_C0_one_over_a2_minus_one_over_a02;
one_over_C0_four_thirds_over_acb_a0cb = 4/3*one_over_C0_one_over_a3_minus_one_over_a03;

Mo_first_term_small = 4*one_over_C0_one_over_a_minus_one_over_a0 + 3*alpha0_sq_small.*one_over_C0_log_a_over_a0;
Mo_second_term_small = 1.5*one_over_C0_alpha0sq_over_asq_a0sq + one_over_C0_four_thirds_over_acb_a0cb;

omega_term_small = 2.5*one_over_C0_alpha0sq_over_asq_a0sq + one_over_C0_four_thirds_over_acb_a0cb;

Mo_first_term = zeros(size(a_0));
Mo_second_term = zeros(size(a_0));
Mo_first_term(find_largeC0) = Mo_first_term_large;
Mo_first_term(find_smallC0) = Mo_first_term_small;

Mo_second_term(find_largeC0) = Mo_second_term_large;
Mo_second_term(find_smallC0) = Mo_second_term_small;

omega_term = zeros(size(a_0));
omega_term(find_largeC0) = omega_term_large;
omega_term(find_smallC0) = omega_term_small;

Mo = 0.125*Mo_first_term + threek2/(16*mu)*(3*c_sq-1).*Mo_second_term + Mo_0;

omega = threek2/(16*mu)*(5*c_sq-1).*omega_term + omega_0;

bigO = -threek2/(8*mu)*c.*omega_term + bigO_0;

figure
plot(C_0,(4./a - 4./a_0),'bo')
hold on
plot(C_0(find_smallC0),4*one_over_C0_one_over_a_minus_one_over_a0.*C0_small,'rx')
set(gca,'Xscale','log')
set(gca,'Yscale','log')

figure
plot(C_0,(log(a./a_0))./C_0,'bo')
hold on
plot(C_0(find_smallC0),one_over_C0_log_a_over_a0,'rx')
set(gca,'Xscale','log')

figure
plot(C_0,log(a./a_0),'ko')
hold on
plot(C_0(find_smallC0),one_over_C0_log_a_over_a0.*C0_small,'rx')
set(gca,'Xscale','log')
% set(gca,'Yscale','log')

figure(11)
hold on
plot(C_0,Mo,'rx')
set(gca,'Xscale','log')

figure(12)
hold on
plot(C_0,omega,'rx')
set(gca,'Xscale','log')

figure(13)
hold on
plot(C_0,bigO,'rx')
set(gca,'Xscale','log')

a_0 = 8e3; %km
e_0 = 0.01;
C_0 = 1e-18;
dt = 1000;
mu = 3.9860e+05;

n_0 = sqrt(mu)*a_0.^(-3/2);
beta_0 = (sqrt(3)/2)*e_0;
beta0_sq = beta_0.^2;

a = a_0/beta0_sq*tan(atan(beta_0)-beta_0*n_0*a_0*C_0*dt)^2;

d1 = -2;
d2 = 1+3*beta0_sq;
d3 = -4/3*beta0_sq*(2+3*beta0_sq);
d4 = 1/3*beta0_sq*(15*beta0_sq^2+15*beta0_sq+2);
d5 = -2/15*beta0_sq*(45*beta0_sq^2+60*beta0_sq+17);

naCdt = n_0*a_0*C_0*dt;

eq35_org = a-a_0
eq35 = a_0*(1+beta0_sq)*(naCdt*d1+naCdt^2*d2+naCdt^3*d3+naCdt^4*d4+naCdt^5*d5)

g1 = -2;
g2 = -1+beta0_sq;
g3 = -2/3*(1+beta0_sq^2);
g4 = 1/6*(1-beta0_sq)*(3*beta0_sq^2+4*beta0_sq+3);
g5 = -2/15*beta0_sq*(3*beta0_sq^4+2*beta0_sq^3+2*beta0_sq+3);

eq36_org = 1/C_0*log(a/a_0)
eq36 = n_0*a_0*(1+beta0_sq)*dt*(naCdt*g1+naCdt^2*g2+naCdt^3*g3+naCdt^4*g4+naCdt^5*g5)

% figure(121)
% hold on
% plot(C_0,mod(Mo,2*pi)-mod(Mo1,2*pi),'rx')
% set(gca,'Xscale','log')
% 
% figure(122)
% hold on
% plot(C_0,omega-omega1,'rx')
% set(gca,'Xscale','log')
% 
% figure(123)
% hold on
% plot(C_0,bigO-bigO1,'rx')
% set(gca,'Xscale','log')

out_oe = [a,e,mod(inc_0, 2*pi),mod(bigO, 2*pi),mod(omega, 2*pi),mod(Mo, 2*pi)];
% find_notreal_inc0 = find(~isreal(inc_0))
% find_notreal_bigO = find(~isreal(bigO))
% find_notreal_omega = find(~isreal(omega))
% find_notreal_Mo = find(~isreal(Mo))
not_real = ~isreal(inc_0) | ~isreal(bigO) | ~isreal(omega) | ~isreal(Mo);
errors(not_real) = 1;
out_oe(not_real,:) = input_oe(not_real,:);
% if isreal(inc_0) && isreal(bigO) && isreal(omega) && isreal(Mo)
%     out_oe=[a;e;mod(inc_0, 2*pi);mod(bigO, 2*pi);mod(omega, 2*pi);mod(Mo, 2*pi)];
% else
%     out_oe=input_oe;
%     error=1;
% end
