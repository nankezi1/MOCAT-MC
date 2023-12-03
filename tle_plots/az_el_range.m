function [az,el,range,x_sight]=az_el_range(r,sight_long,sight_lat,sight_height,gst);

%function [az,el,range]=az_el_range(r,sight_long,sight_lat,sight_height,gst);
%
% This function computes the azimuth, elevation and range to a satellite.
%
%  The inputs are:
%               r = position of satellite in ECI coordinates (mx3)
%      sight_long = longitude of sensor sight in degrees (1x1)
%       sight_lat = latitude of sensor sight in degrees (1x1)
%    sight_height = latitude of sensor sight in km (1x1)
%             gst = Greenwich Mean Sidereal Time (deg) (mx1)
% 
%  The outputs are:
%       az = azimuth in degrees (mx1), note from 0 to 360 degrees
%       el = elevation in degrees (mx1), note from -90 to 90 degrees
%    range = range in km (mx1)

% John L. Crassidis 5/9/11

% Convert to Radians
gst=gst*pi/180;
sight_long=sight_long*pi/180;
sight_lat=sight_lat*pi/180;

% Earth Parameters
rearth=6378.1363;ecc_earth2=0.006694385;

% Flat Earth Model from Vallado (third edition), p. 149
cos_comp=rearth/sqrt(1-ecc_earth2*sin(sight_lat)^2);
sin_comp=rearth*(1-ecc_earth2)/sqrt(1-ecc_earth2*sin(sight_lat)^2);
r_delta=(cos_comp+sight_height)*cos(sight_lat);
r_k=(sin_comp+sight_height)*sin(sight_lat);
r_sight=(r_delta^2+r_k^2)^(0.5);

% Sideral Time of Sensor Sight
theta=gst+sight_long;

% Vector from Sensor Sight to Satellite in Inertial Coordinates
rho1_sight=r(:,1)-r_sight*cos(sight_lat).*cos(theta);
rho2_sight=r(:,2)-r_sight*cos(sight_lat).*sin(theta);
rho3_sight=r(:,3)-r_sight*sin(sight_lat);
rho_sight=[rho1_sight rho2_sight rho3_sight];
range=(rho1_sight.^2+rho2_sight.^2+rho3_sight.^2).^(0.5);

% Vector from Sensor Sight to Satellite in Up-East-North Coordinates
rho_u_sight=cos(sight_lat)*cos(theta).*rho1_sight+cos(sight_lat)*sin(theta).*rho2_sight+sin(sight_lat)*rho3_sight;
rho_e_sight=-sin(theta).*rho1_sight+cos(theta).*rho2_sight;
rho_n_sight=-sin(sight_lat)*cos(theta).*rho1_sight-sin(sight_lat)*sin(theta).*rho2_sight+cos(sight_lat)*rho3_sight;

% Azimuth and Elevation of sight
az=rem(atan2(rho_e_sight,rho_n_sight)*180/pi+360,360);
el=asin(rho_u_sight./range)*180/pi;
x_sight=[rho_n_sight,rho_e_sight,rho_u_sight];
