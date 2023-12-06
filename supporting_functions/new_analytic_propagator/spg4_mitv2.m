function [sat,r_pos,v_vel]=spg4_mitv2(sat,time)
%constants
CK2 = 5.413080e-4;
CK4 = 0.62098875e-6;
E6A = 1.0e-6;

%QOMS2T = 1.88027916e-9;
 QOMS2T =     1.880276800610897e-09;

%S = 1.01222928;
S=1.012229276354522e+00;

%XJ3 = -0.253881e-5;
XJ3 =-2.532153060000000e-06;

%XKE = 0.743669161e-1;
XKE =7.436685316871385e-02;

XKMPER = 6378.137;
XMNPDA = 1440.0;
AE = 1.0;
% earth flattening
F = 1/298.257223563;



% self.tle["inclination"] = np.deg2rad(self.tle["inclination"])
% self.tle["right_ascension"] = np.deg2rad(self.tle["right_ascension"])
% self.tle["arg_perigee"] = np.deg2rad(self.tle["arg_perigee"])
% self.tle["mean_anomaly"] = np.deg2rad(self.tle["mean_anomaly"])

% self.tle["mean_motion"] *= (np.pi * 2 / XMNPDA)
% self.tle["mean_motion_derivative"] *= np.pi * 2 / XMNPDA ^ 2
% self.tle["mean_motion_sec_derivative"] *= np.pi * 2 / XMNPDA ^ 3
% self.tle["bstar"] *= AE

% n_0 = sat.no;%self.tle["mean_motion"]
% k_e = XKE;
% k_2 = CK2;
% i_0 =sat.inclo; %self.tle["inclination"]
% e_0 =sat.ecco; %self.tle["excentricity"]
% 
% a_1 = (k_e / n_0) ^ (2.0/3);
% delta_1 = ((3/2.0) * (k_2 / a_1^2) * ((3 * cos(i_0)^2 - 1) /...
%                                       (1 - e_0^2)^(2.0/3)));
% 
% a_0 = a_1 * (1 - delta_1/3 - delta_1^2 - (134.0/81) * delta_1^3);
% 
% delta_0 = ((3/2.0) * (k_2 / a_0^2) * ((3 * np.cos(i_0)^2 - 1) /...
%                                       (1 - e_0^2)^(2.0/3)));
% 
% % original mean motion
% n_0pp = n_0 / (1 + delta_0);
%self.tle["original_mean_motion"] = n_0pp

% % semi major axis
% a_0pp = a_0 / (1 - delta_0);
% %self.tle["semi_major_axis"] = a_0pp
% 
% %self.tle["period"] = np.pi * 2 / n_0pp
% 
% perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER;%self.tle["perigee"]

%now = self.tle["epoch"]

% self.tle["right_ascension_lon"] = (self.tle["right_ascension"]
%                                    - gmst(now))
% 
% if self.tle["right_ascension_lon"] > np.pi:
%     self.tle["right_ascension_lon"] -= 2 * np.pi


% no=sat.no;
% Omegao=sat.nodeo;
% eo=sat.ecco;
% io=sat.inclo;
% Mo=sat.mo;
% omegao=sat.argpo;
% Bstar=sat.bstar;

%perigee = sat.a*(1-sat.ecco);%["perigee"]
a_0pp = sat.a_0pp;%["semi_major_axis"]
e_0 = sat.ecco;%["excentricity"]
i_0 = sat.inclo;%["inclination"]
n_0pp = sat.n_0pp;%["original_mean_motion"]
k_2 = CK2;
k_4 = CK4;
k_e = XKE;
bstar = sat.bstar;%["bstar"]
w_0 = sat.argpo;%["arg_perigee"]
M_0 = sat.mo;%["mean_anomaly"]
W_0 = sat.nodeo;%["right_ascension"]
t_0 = 0;
A30 = -XJ3 * AE^3;

perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER;

if perigee < 98
    s = 20/XKMPER + AE;
    qoms2t = (QOMS2T^0.25 + S - s)^4;
elseif perigee < 156
    s = a_0pp * (1 - e_0) - S + AE ;
    qoms2t = (QOMS2T ^ 0.25 + S - s)^4;
else
    qoms2t = QOMS2T;
    s = S;
end

theta = cos(i_0);
xi = 1 / (a_0pp - s);
beta_0 = sqrt(1 - e_0 ^ 2);
eta = a_0pp * e_0 * xi;

C_2 = (qoms2t * xi^4 * n_0pp * (1 - eta^2)^(-3.5) *...
       (a_0pp * (1 + 1.5 * eta^2 + 4 * e_0 * eta + e_0 * eta^3) +...
        1.5 * (k_2 * xi) / (1 - eta^2) * (-0.5 + 1.5 * theta^2)*...
        (8 + 24 * eta^2 + 3 * eta^4)));

C_1 = bstar * C_2;

C_3 = (qoms2t * xi ^ 5 * A30 * n_0pp * AE * sin(i_0) / (k_2 * e_0));

coef = 2 * qoms2t * xi^4 * a_0pp * beta_0^2*(1-eta^2)^(-7/2.0);

C_4 = (coef * n_0pp *...
       ((2 * eta * (1 + e_0 * eta) + e_0/2.0 + (eta^3)/2.0) -...
        2 * k_2 * xi / (a_0pp * (1 - eta^2)) *...
        (3*(1-3*theta^2) *...
         (1 + (3*eta^2)/2.0 - 2*e_0*eta - e_0*eta^3/2.0) +...
         3/4.0*(1-theta^2)*...
         (2*eta^2 - e_0*eta - e_0*eta^3)*cos(2*w_0))));

C_5 = coef * (1 + 11/4.0 * eta * (eta + e_0) + e_0 * eta^3);
D_2 = 4 * a_0pp * xi * C_1^2;
D_3 = 4/3.0 * a_0pp * xi^2 * (17*a_0pp + s) * C_1^3;
D_4 = 2/3.0 * a_0pp * xi^3 * (221*a_0pp + 31*s) * C_1^4;

% Secular effects of atmospheric drag and gravitation
%dt = _days(current_time - t_0) * XMNPDA
dt=time-t_0;

M_df = (M_0 + (1 +...
               3*k_2*(-1 + 3*theta^2)/(2*a_0pp^2 * beta_0^3) +...
               3*k_2^2*(13 - 78*theta^2 + 137*theta^4)/...
               (16*a_0pp^4*beta_0^7))*...
        n_0pp*dt);
w_df = (w_0 + (-3*k_2*(1 - 5*theta^2)/(2*a_0pp^2*beta_0^4) +...
               3 * k_2^2 * (7 - 114*theta^2 + 395*theta^4)/...
               (16*a_0pp*beta_0^8) +...
               5*k_4*(3-36*theta^2+49*theta^4)/...
               (4*a_0pp^4*beta_0^8))*...
        n_0pp*dt);
W_df = (W_0 + (-3*k_2*theta/(a_0pp^2*beta_0^4) +...
               3*k_2^2*(4*theta- 19*theta^3)/(2*a_0pp^4*beta_0^8) +...
               5*k_4*theta*(3-7*theta^2)/(2*a_0pp^4*beta_0^8))*...
        n_0pp*dt);
deltaw = bstar * C_3 * cos(w_0)*dt;
deltaM = (-2/3.0 * qoms2t * bstar * xi^4 * AE / (e_0*eta) *...
          ((1 + eta * cos(M_df))^3 - (1 + eta * cos(M_0))^3));
M_p = M_df + deltaw + deltaM;
w = w_df - deltaw - deltaM;
W = (W_df - 21/2.0 * (n_0pp * k_2 * theta)/(a_0pp^2 * beta_0^2) *...
     C_1 * dt^2);

e = (e_0 -...
     bstar * C_4 * dt -...
     bstar * C_5 * (sin(M_p) - sin(M_0)));

a = a_0pp * (1 - C_1 * dt - D_2 * dt^2 - D_3 * dt^3 - D_4 * dt^4)^2;
L = M_p + w + W + n_0pp * (3/2.0 * C_1 * dt^2 +...
                           (D_2 + 2 * C_1 ^ 2) * dt^3 +...
                           1/4.0 *...
                           (3*D_3 + 12*C_1*D_2 + 10*C_1^3)*dt^4 +...
                           1.0/5 * (3*D_4 + 12*C_1*D_3 + 6*D_2^2 +...
                                    30*C_1^2*D_2 + 15*C_1^4)*dt^5);
beta = sqrt(1 - e^2);
n = k_e / (a ^ (3/2.0));

% Long-period periodic terms
a_xN = e * cos(w);
a_yNL = A30 * sin(i_0) / (4.0 * k_2 * a * beta^2);
L_L = a_yNL/2 * a_xN * ((3 + 5 * theta) / (1 + theta));
L_T = L + L_L;
a_yN = e * sin(w) + a_yNL;

U = mod((L_T - W),pi*2); % (np.pi * 2)

Epw = U;
for k =1:100
    DeltaEpw = ((U - a_yN * cos(Epw) + a_xN  * sin(Epw) - Epw) /...
                (-a_yN * sin(Epw) - a_xN * cos(Epw) + 1));
    Epw = Epw + DeltaEpw;
    if DeltaEpw < 10e-14
        break
    end
end
% preliminary quantities for short-period periodics

ecosE = a_xN * cos(Epw) + a_yN * sin(Epw);
esinE = a_xN * sin(Epw) - a_yN * cos(Epw);

e_L = (a_xN^2 + a_yN^2)^(0.5);
p_L = a * (1 - e_L^2);
r = a * (1 - ecosE);
rdot = k_e * sqrt(a)/r * esinE;
rfdot = k_e * sqrt(p_L) / r;
cosu = a / r * (cos(Epw) - a_xN +...
                (a_yN * (esinE) / (1 + sqrt(1 - e_L^2))));
sinu = a / r * (sin(Epw) - a_yN +...
                (a_xN * (esinE) / (1 + sqrt(1 - e_L^2))));
u = atan2(sinu, cosu);


cos2u = cos(2*u);
sin2u = sin(2*u);

Deltar = k_2/(2*p_L) * (1 - theta^2) * cos2u;
Deltau = -k_2/(4*p_L^2) * (7*theta^2 - 1) * sin2u;
DeltaW = 3*k_2 * theta / (2 * p_L^2) * sin2u;
Deltai = 3*k_2 * theta / (2 * p_L^2) * cos2u * sin(i_0);
Deltardot = - k_2 * n / p_L * (1 - theta^2) * sin2u;
Deltarfdot = k_2 * n / p_L * ((1 - theta^2) * cos2u -...
                              3/2.0 * (1 - 3*theta^2));

% osculating quantities

r_k = r * (1 - 3/2.0 * k_2 * sqrt(1 - e_L^2)/p_L^2 *...
           (3 * theta^2 - 1)) + Deltar;
u_k = u + Deltau;
W_k = W + DeltaW;
i_k = i_0 + Deltai;
rdot_k = rdot + Deltardot;
rfdot_k = rfdot + Deltarfdot;

M_x = -sin(W_k) * cos(i_k);
M_y = cos(W_k) * cos(i_k);
M_z = sin(i_k);

N_x = cos(W_k);
N_y = sin(W_k);
N_z = 0;

U_x = M_x * sin(u_k) + N_x * cos(u_k);
U_y = M_y * sin(u_k) + N_y * cos(u_k);
U_z = M_z * sin(u_k) + N_z * cos(u_k);

V_x = M_x * cos(u_k) - N_x * sin(u_k);
V_y = M_y * cos(u_k) - N_y * sin(u_k);
V_z = M_z * cos(u_k) - N_z * sin(u_k);


r_x = r_k * U_x;
r_y = r_k * U_y;
r_z = r_k * U_z;

rdot_x = rdot_k * U_x + rfdot_k * V_x;
rdot_y = rdot_k * U_y + rfdot_k * V_y;
rdot_z = rdot_k * U_z + rfdot_k * V_z;

r_pos=XKMPER*[r_x, r_y, r_z]'; 
v_vel=[rdot_x, rdot_y, rdot_z]';

% no=sat.no;
% Omegao=sat.nodeo;
% eo=sat.ecco;
% io=sat.inclo;
% Mo=sat.mo;
% omegao=sat.argpo;
% Bstar=sat.bstar;

sat.no=n;
sat.nodeo=W;
sat.ecco=e;
sat.inclo=i_0;
sat.argpo=w;
sat.mo=M_p;
sat.n_0pp=n;
sat.a_0pp=a;

