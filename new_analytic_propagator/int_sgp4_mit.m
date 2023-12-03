function [sat]=int_sgp4_mit(sat)

%constants
CK2 = 5.413080e-4;
CK4 = 0.62098875e-6;
E6A = 1.0e-6;
QOMS2T = 1.88027916e-9;
S = 1.01222928;
XJ3 = -0.253881e-5;
XKE = 0.743669161e-1;
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

n_0 = sat.no;%self.tle["mean_motion"]
k_e = XKE;
k_2 = CK2;
i_0 =sat.inclo; %self.tle["inclination"]
e_0 =sat.ecco; %self.tle["excentricity"]

a_1 = (k_e / n_0) ^ (2.0/3);
delta_1 = ((3/2.0) * (k_2 / a_1^2) * ((3 * cos(i_0)^2 - 1) /...
                                      (1 - e_0^2)^(2.0/3)));

a_0 = a_1 * (1 - delta_1/3 - delta_1^2 - (134.0/81) * delta_1^3);

delta_0 = ((3/2.0) * (k_2 / a_0^2) * ((3 * cos(i_0)^2 - 1) /...
                                      (1 - e_0^2)^(2.0/3)));

% original mean motion
n_0pp = n_0 / (1 + delta_0);
%self.tle["original_mean_motion"] = n_0pp

% semi major axis
a_0pp = a_0 / (1 - delta_0);
%self.tle["semi_major_axis"] = a_0pp

%self.tle["period"] = np.pi * 2 / n_0pp

perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER;%self.tle["perigee"]

%now = self.tle["epoch"]

% self.tle["right_ascension_lon"] = (self.tle["right_ascension"]
%                                    - gmst(now))
% 
% if self.tle["right_ascension_lon"] > np.pi:
%     self.tle["right_ascension_lon"] -= 2 * np.pi

sat.n_0pp=n_0pp;
sat.a_0pp=a_0pp;