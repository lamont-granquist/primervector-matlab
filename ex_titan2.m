close all; clear all; clc;
format longG;

%
% BODY CONSTANTS
%

body.mu = 3.986004418e+14;         % m^3/s^2; earth
body.r  = 6.371e+6;                % m; earth radius
body.w  = 7.2921159e-5 * [0,0,1]'; % rad/s; sidereal angular

%
% VEHICLE
%

payload = 3580;   % kg

stages(1).m0     = 149600 + payload;  % kg
stages(1).thrust = 2 * 1097.2 * 1000; % N
stages(1).isp    = 296;               % sec
stages(1).bt     = 156;               % sec

stages(2).m0     = 28400 + payload;   % kg
stages(2).thrust = 443.7 * 1000;        % N
stages(2).isp    = 315;               % sec
%stages(2).bt     = 180;               % sec
stages(2).bt = 166.753090665355;

%
% BOUNDARY CONDITIONS
%

incT = deg2rad(28.608);
PeA = 185e+3;
ApA = 185e+3;
gammaT = deg2rad(0);

bcfun = @(xf) BCflightangle4constraintPeAApA(xf, body, PeA, ApA, incT, gammaT);

% launchsite latitude
lat = deg2rad(28.608);
lng = 0;

[rf vf] = launch(stages, body, bcfun, lat, lng, incT)

%
% SOME OUTPUT
%

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf', vf', body.mu);

orbit.semi_major_axis = a / 1000;
orbit.eccentricity = eMag;
orbit.PeR = (1 - eMag) * a;
orbit.ApR = (1 + eMag) * a;
orbit.PeA = orbit.PeR - body.r;
orbit.ApA = orbit.ApR - body.r;
orbit.inclination = rad2deg(i);
orbit.LAN = rad2deg(O);
orbit.argument_of_periapsis = rad2deg(o);
orbit.true_anomaly = rad2deg(nu);
orbit.semi_latus_rectum = p / 1000;

disp(orbit);

