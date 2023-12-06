function [out_oe,error]=analytic_propagation(input_oe,param)

error=0;
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

if strcmpi(param.density_profile,'JB2008')

    % if statement ensure we get reasonable values for rho. Rho=0 break the
    % model below. Also altitude greater than 2300 km breaks the model.
    if (input_oe(1,1)-re)<2000 && (input_oe(1,1)-re)>200% changes from 200-1100
    
        %rho_0 = interp2(param.dens_times,param.alt,param.dens_value,param.jd,input_oe(1,1)-re)*(1000)^3;
    %     rho_0 =qinterp2((param.alt)',(param.dens_times)',param.dens_value',input_oe(1,1)-re, param.jd,1)*(1000)^3;
        %     F107=param.F107;
        %     Ap=param.Ap;
        %     % density model
        %     h=input_oe(1,1);
        %     T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
        %     m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
        %     H = T / m;
        %     rho_0 = 6e-10* exp ( - ( input_oe(1,1)-re- 175 ) / H ) *(1000)^3;%units kg/km^3
    %     rho_old = qinterp2((param.alt)',(param.dens_times)',param.dens_value',input_oe(1,1)-re, param.jd,1)*(1000)^3;
    %     rho_grid =param.dens_grid(input_oe(1,1)-re, param.jd);
        rho_0 = lininterp2(param.alt(:,1),param.dens_times(1,:),param.dens_value,input_oe(1,1)-re, param.jd)*(1000)^3;
        % https://www.mathworks.com/matlabcentral/fileexchange/28376-faster-linear-interpolation
    
    %     fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
    else
        if (input_oe(1,1)-re)>1100
            %         aa=1100+re;
            %         F107=param.F107;
            %         Ap=param.Ap;
            %         % density model
            %         h=aa;
            %         T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
            %         m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
            %         H = T / m;
            %         rho_0 = 6e-10* exp ( - ( aa-re- 175 ) / H ) *(1000)^3;%units kg/km^3
            % rho_0 = interp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
            rho_0 = lininterp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
    %         rho_old = interp1(param.dens_times(1,:),param.dens_value(end,:),param.jd)*(1000)^3;
    %         fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
    
        else
            %         aa=200+re;
            %         F107=param.F107;
            %         Ap=param.Ap;
            %         % density model
            %         h=aa;
            %         T = 900 + 2.5 *( F107 - 70 ) + 1.5* Ap;
            %         m = 27 - 0.012* ( h-re - 200 ); %180 < h(km) < 500
            %         H = T / m;
            %         rho_0 = 6e-10* exp ( - ( aa-re- 175 ) / H ) *(1000)^3;%units kg/km^3
            % rho_0 = interp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
            rho_0 = lininterp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
    %         rho_old = interp1(param.dens_times(1,:),param.dens_value(1,:),param.jd)*(1000)^3;
    %         fprintf('old vs new: %0.3e \t %0.3e \n', rho_0, rho_old);
        end
    end

elseif strcmpi(param.density_profile,'static')
    rho_0 = densityexp(input_oe(1,1)-re)*(1000)^3;
end
%constants
%C_0=param.C_0*rho_0;%1/2*Cd*s_ref/m*rho_0

%check these unit unit of 1/re and using 1/(6378.137 km) to covert to km
%units

% rho_reference for Bstar is 0.157 in units of kg/(m^2 * re).
% Using (1000^2 m^2/1 km^2) to convert A from m^2 to km^2, and dividing by 0.157 to scale
% according to newly computed density rho_0 (kg/km^3)
C_0 = max((abs(param.Bstar)/(1e6*0.157))*rho_0,1e-16); 
% C_0=(abs(param.Bstar)/((0.157/re)*1000^2))*rho_0;

%C_0=param.C0*rho_0;
k2=mu*J2*re^2/2;

t=param.t;
t_0=param.t_0;

%initial conditions
a_0=input_oe(1,1);
e_0=input_oe(2,1);
inc_0=input_oe(3,1);
bigO_0=input_oe(4,1);
omega_0=input_oe(5,1);
Mo_0=input_oe(6,1);

c=cos(inc_0);
n_0=sqrt(mu)*a_0^(-3/2);

% p_0=a_0*(1-e_0^2);

% Cu=(3/4)*n_0*J2*(re/p_0)^2;

alpha_0=e_0/sqrt(a_0);

% if e_0<100000e-2

beta_0=(sqrt(3)/2)*e_0;

tan_coeff = max((tan(atan(beta_0)-beta_0*n_0*a_0*C_0*(t-t_0))),0);
a=(a_0/beta_0^2)*tan_coeff^2;

e=(1/(sqrt(3)/2))*tan_coeff;

Mo= (1/8)*(1/C_0)*(4/a+3*alpha_0^2*log(a/a_0))+...
    -(1/8)*(1/C_0)*(4/a_0+3*alpha_0^2*log(a_0/a_0))+...
    (3*k2*(3*c^2-1))/(16*mu)*(1/C_0)*((3*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    -(3*k2*(3*c^2-1))/(16*mu)*(1/C_0)*((3*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+...
    Mo_0;

omega=(3*k2*(5*c^2-1))/(16*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    -(3*k2*(5*c^2-1))/(16*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+omega_0;

bigO=-(3*k2*c)/(8*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a^2+4/(3*a^3))+...
    (3*k2*c)/(8*mu)*(1/C_0)*((5*alpha_0^2/2)*1/a_0^2+4/(3*a_0^3))+bigO_0;

% else
%     %out_oe_rk=rk_mean_elements(input_oe,C_0,t-t_0);
%
%     options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%     [~,y] = ode45(@mean_elements_ode,[t_0 t],input_oe,options,C_0);
%     out_oe_rk=y(end,:)';
%
%     a=out_oe_rk(1,1);
%     e=out_oe_rk(2,1);
%     inc=out_oe_rk(3,1);
%     bigO=out_oe_rk(4,1);
%     omega=out_oe_rk(5,1);
%     Mo=out_oe_rk(6,1);
%  end

if isreal(inc_0) && isreal(bigO) && isreal(omega) && isreal(Mo)
    out_oe=[a;e;mod(inc_0, 2*pi);mod(bigO, 2*pi);mod(omega, 2*pi);mod(Mo, 2*pi)];
else
    out_oe=input_oe;
    error=1;
end
