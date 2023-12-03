function varargout = sun_position(JD, rsat)

% Algorithm 29 in Vallado
 
% Julian Centuries 
T_UT1 = (JD - 2451545.0)/36525;

% Mean Longitude of the Sun
meanlong_deg = 280.460 + 36000.77*T_UT1;
meanlong_deg = rem(meanlong_deg, 360);

if meanlong_deg < 0
    meanlong_deg = meanlong_deg + 360;
end

% Mean anomoly of the sun
meananomaly_deg = 357.5277233 + 35999.05034*T_UT1;
meananomaly_deg = rem(meananomaly_deg, 360);

if meananomaly_deg < 0
    meananomaly_deg = meananomaly_deg + 360;
end

% Ecliptic Parameters
lambda_ecliptic_deg = meanlong_deg + 1.914666471*sin(meananomaly_deg*pi/180) + 0.019994643*sin(2*meananomaly_deg*pi/180);
obliquity_deg = 23.439291 - 0.0130042*T_UT1;

% Sun Distance in au
r_sun = 1.000140612 - 0.016708617*cos(meananomaly_deg*pi/180) - 0.000139589*cos(2*meananomaly_deg*pi/180);

% Position Vector from the Earth to the Sun in ECI (au)
rvec_sun = r_sun*[cos(lambda_ecliptic_deg*pi/180) 
                  cos(obliquity_deg*pi/180)*sin(lambda_ecliptic_deg*pi/180) 
                  sin(obliquity_deg*pi/180)*sin(lambda_ecliptic_deg*pi/180)];
              
% Sun Right Ascension and Declanation
delta_sun = asin(rvec_sun(3)/r_sun)*180/pi;
alpha_sun = atan2(rvec_sun(2)/r_sun, rvec_sun(1)/r_sun)*180/pi;

if alpha_sun < 0
    alpha_sun = alpha_sun + 360;
end

delta_sun = delta_sun*pi/180;
alpha_sun = alpha_sun*pi/180;

%% Determine if the satellite is shadowed

% Algorithm 35 in Vallado
% Convert Position vector of sun and satellite to Earth Radii
rsun = rvec_sun * 149597870.0/6378.137;      % Converts au to ER
rsun_out = rsun;
rsat = rsat * 1/6378.137;                    % Converts km to ER

% Unit vector from Sat to Sun
uvec_satsun = (rsun - rsat)/norm(rsun - rsat);

% Distance From Satellite to Sun
d_satsun = norm(rsun - rsat) * 6378.137/149597870.0;     % Expressed in au

% Scale the z-components so to account for shadowing of a non-spherical planet
e2 = 0.006694385;
rsun(3,1) = rsun(3,1)/sqrt(1-e2);
rsat(3,1) = rsat(3,1)/sqrt(1-e2);

% Parametric Value
tau_min = (rsat'*rsat - rsat'*rsun)/(rsat'*rsat + rsun'*rsun - 2*rsat'*rsun); 

if tau_min < 0 || tau_min > 1
    in_shadow = 0;
else
    c2 = (1-tau_min)*(rsat'*rsat) + (rsat'*rsun)*tau_min;
    if c2 >= 1 % Not in the shadow
        in_shadow = 0;
    else % In the shadow
        in_shadow = 1;
    end
end

varargout{1} = d_satsun;
varargout{2} = uvec_satsun;
varargout{3} = in_shadow;
varargout{4} = rsun_out;
varargout{5} = alpha_sun;
varargout{6} = delta_sun;
end