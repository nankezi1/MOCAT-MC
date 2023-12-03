function [out_oe,errors]=analytic_propagation_vec_ecc(input_oe,param)

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
a_0 = input_oe(:,1);
a_minus_re = a_0-re;

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

% C_0 = max((param.Bstar/(1e6*0.157)).*rho_0,1e-16);
C_0 = logspace(-8,-19,numel(a_0))';

%C_0=param.C0*rho_0;
k2_over_mu = J2*re^2/2;

t = param.t;
t_0 = param.t_0;

%initial conditions
e_0 = input_oe(:,2);
inc_0 = input_oe(:,3);
bigO_0 = input_oe(:,4);
omega_0 = input_oe(:,5);
Mo_0 = input_oe(:,6);

figure(10)
plot(1:numel(a_0),sort(e_0),'kx')

e_0 = logspace(log10(max(e_0)),-40,numel(a_0))';

c = cos(inc_0);
c_sq = c.^2;

tic
%%%%%No eccentricity exception
n_0 = sqrt(mu)*a_0.^(-3/2);

% p_0=a_0.*(1-e_0.^2);

% Cu=(3/4)*n_0*J2.*(re./p_0).^2;

alpha_0 = e_0./sqrt(a_0);
alpha0_sq = alpha_0.^2;

beta_0 = (sqrt(3)/2)*e_0;
beta0_sq = beta_0.^2;

tan_atan_beta0 = max(tan(atan(beta_0)-beta_0.*n_0.*a_0.*C_0*(t-t_0)),0); %place lower limit on eccentricity and semi-major axis reduction
a = (a_0./beta0_sq).*tan_atan_beta0.^2;
e = (2/sqrt(3))*tan_atan_beta0;

%Save some variables to avoid repetition of operations
a_sq = a.^2;
four_thirds_over_a_cb = 4/3./(a_sq.*a);
a0_sq = a_0.^2;
four_thirds_over_a0_cb = 4/3./(a0_sq.*a_0);
alpha0sq_over_asq = alpha0_sq./a_sq;
alpha0sq_over_a0sq = alpha0_sq./a0_sq;

Mo = 0.125*(4./a - 4./a_0 + 3*alpha0_sq.*log(a./a_0))./C_0 + 3*k2_over_mu/16*(3*c_sq-1).*(1.5*(alpha0sq_over_asq-alpha0sq_over_a0sq)+four_thirds_over_a_cb-four_thirds_over_a0_cb)./C_0 + Mo_0;

five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 = (2.5*(alpha0sq_over_asq-alpha0sq_over_a0sq)+four_thirds_over_a_cb-four_thirds_over_a0_cb)./C_0;

omega = 3*k2_over_mu/16*(5*c_sq-1).*five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 + omega_0;

bigO = -3*k2_over_mu/8*c.*five_a0sq_over2_tau2_plus_4thirds_over_tau3_overC0 + bigO_0;
toc

figure(1)
plot(e_0,Mo,'bo')
set(gca,'XScale','log')

figure(2)
plot(e_0,omega,'bo')
set(gca,'XScale','log')

figure(3)
plot(e_0,bigO,'bo')
set(gca,'XScale','log')

tic
%%%%%%%%%
check_e0 = e_0<1e-3;

%Large eccentricity
e0_elarge = e_0(~check_e0);
a0_elarge = a_0(~check_e0);

n0_elarge = sqrt(mu)*a0_elarge.^(-3/2);

alpha0_elarge = e0_elarge./sqrt(a0_elarge);
alpha0_sq_elarge = alpha0_elarge.^2;

beta_0 = (sqrt(3)/2)*e0_elarge;
beta0_sq = beta_0.^2;

C0_elarge = C_0(~check_e0);

tan_atan_beta0 = max(tan(atan(beta_0)-beta_0.*n0_elarge.*a0_elarge.*C0_elarge*(t-t_0)),0); %place lower limit on eccentricity and semi-major axis reduction
a_elarge = (a0_elarge./beta0_sq).*tan_atan_beta0.^2;
e_elarge = (2/sqrt(3))*tan_atan_beta0;

%Save some variables to avoid repetition of operations
a_sq = a_elarge.^2;
four_thirds_over_a_cb = 4/3./(a_sq.*a_elarge);
a0_sq = a0_elarge.^2;
four_thirds_over_a0_cb = 4/3./(a0_sq.*a0_elarge);
alpha0sq_over_asq = alpha0_sq_elarge./a_sq;
alpha0sq_over_a0sq = alpha0_sq_elarge./a0_sq;

Mo_first_term_elarge = 0.125*(4./a_elarge - 4./a0_elarge + 3*alpha0_sq_elarge.*log(a_elarge./a0_elarge))./C0_elarge;
Mo_second_term_elarge = 3/16*(1.5*(alpha0sq_over_asq-alpha0sq_over_a0sq)+four_thirds_over_a_cb-four_thirds_over_a0_cb)./C0_elarge;

omega_term_elarge = 3/16*(2.5*(alpha0sq_over_asq-alpha0sq_over_a0sq)+four_thirds_over_a_cb-four_thirds_over_a0_cb)./C0_elarge;

%Small eccentricity
a0_esmall = a_0(check_e0);
e0_esmall = e_0(check_e0);

n0_esmall = sqrt(mu)*a0_esmall.^(-3/2);

C0_esmall = C_0(check_e0);

one_minus_Cnadt = max((1-C0_esmall.*n0_esmall.*a0_esmall*(t-t_0)),0);
a_esmall = a0_esmall.*one_minus_Cnadt.^2;
e_esmall = e0_esmall.*one_minus_Cnadt;

one_over_C0_one_over_a_cube_C0_large = (1./a_esmall.^3 - 1./a0_esmall.^3)./C0_esmall;

Mo_first_term_esmall = 0.5*(1./a_esmall - 1./a0_esmall)./C0_esmall;
Mo_second_term_esmall = 0.25*one_over_C0_one_over_a_cube_C0_large;

omega_term_esmall = Mo_second_term_esmall;

%Combine vectors
a = zeros(size(a_0));
a(~check_e0) = a_elarge;
a(check_e0) = a_esmall;

e = zeros(size(a_0));
e(~check_e0) = e_elarge;
e(check_e0) = e_esmall;

Mo_first_term = zeros(size(a_0));
Mo_second_term = zeros(size(a_0));
Mo_first_term(~check_e0) = Mo_first_term_elarge;
Mo_first_term(check_e0) = Mo_first_term_esmall;

Mo_second_term(~check_e0) = Mo_second_term_elarge;
Mo_second_term(check_e0) = Mo_second_term_esmall;

omega_term = zeros(size(a_0));
omega_term(~check_e0) = omega_term_elarge;
omega_term(check_e0) = omega_term_esmall;

Mo = Mo_first_term + k2_over_mu*(3*c_sq-1).*Mo_second_term + Mo_0;

omega = k2_over_mu*(5*c_sq-1).*omega_term + omega_0;

bigO = -2*k2_over_mu*c.*omega_term + bigO_0;
toc

figure(1)
hold on
plot(e_0,Mo,'rx')

figure(2)
hold on
plot(e_0,omega,'rx')

figure(3)
hold on
plot(e_0,bigO,'rx')

% c=cos(inc_0);
% n_0=sqrt(mu)*a_0.^(-3/2);
% 
% % p_0=a_0.*(1-e_0.^2);
% 
% % Cu=(3/4)*n_0*J2.*(re./p_0).^2;
% 
% alpha_0=e_0./sqrt(a_0);
% 
% % if e_0<100000e-2
% 
% beta_0=(sqrt(3)/2)*e_0;
% 
% a=(a_0./beta_0.^2).*(tan(atan(beta_0)-beta_0.*n_0.*a_0.*C_0*(t-t_0))).^2;
% 
% e=(1/(sqrt(3)/2)).*tan(atan(beta_0)-beta_0.*n_0.*a_0.*C_0*(t-t_0));
% 
% Mo= (1/8).*(1./C_0).*(4./a+3.*alpha_0.^2.*log(a./a_0))+...
%     -(1/8)*(1./C_0).*(4./a_0+3*alpha_0.^2.*log(a_0./a_0))+...
%     (3*k2*(3*c.^2-1))/(16*mu).*(1./C_0).*((3*alpha_0.^2/2)*1./a.^2+4./(3*a.^3))+...
%     -(3*k2*(3*c.^2-1))/(16*mu).*(1./C_0).*((3*alpha_0.^2/2)*1./a_0.^2+4./(3*a_0.^3))+...
%     Mo_0;
% 
% omega=(3*k2*(5*c.^2-1))/(16*mu).*(1./C_0).*((5*alpha_0.^2/2)*1./a.^2+4./(3*a.^3))+...
%     -(3*k2*(5*c.^2-1))/(16*mu).*(1./C_0).*((5*alpha_0.^2/2)*1./a_0.^2+4./(3*a_0.^3))+omega_0;
% 
% bigO=-(3*k2*c)/(8*mu).*(1./C_0).*((5*alpha_0.^2/2)*1./a.^2+4./(3*a.^3))+...
%     (3*k2*c)/(8*mu).*(1./C_0).*((5*alpha_0.^2/2)*1./a_0.^2+4./(3*a_0.^3))+bigO_0;

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
