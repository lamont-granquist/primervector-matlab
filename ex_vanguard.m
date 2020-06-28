close all; clear classes; clear all; clc;
format longG;
format compact;

%
% BODY CONSTANTS
%

body.mu = 3.986004418e+14;         % m^3/s^2; earth
body.r  = 6.371e+6;                % m; earth radius
body.w  = 7.2921159e-5 * [0,0,1]'; % rad/s; sidereal angular

%
% VEHICLE
%

% Vanguard 3
payload = 23.7;   % kg

% SLV-7
stages(1).m0       = 8090 + 1990 + 194 + payload; % kg
stages(1).thrust   = 125 * 1000;                  % N
stages(1).isp      = 248;                         % sec
stages(1).bt       = 144;                         % sec

% AJ10-37
stages(2).m0       = 1990 + 194 + payload; % kg
stages(2).thrust   = 32.6 * 1000;          % N
stages(2).isp      = 261;                  % sec
stages(2).bt       = 120;                  % sec
stages(2).bt_free  = true;

% X248-A2
stages(3).m0       = 194 + payload; % kg
stages(3).thrust   = 11.6 * 1000;   % N
stages(3).isp      = 230;           % sec
stages(3).bt       = 31;            % sed
stages(3).unguided = true;

%
% BOUNDARY CONDITIONS
%

incT = deg2rad(33.4);
PeA = 512e+3;
ApA = 3743e+3;
%AttR = 185e+3;
%LANT = deg2rad(270);
%ArgPT = deg2rad(135);

%bcfun = @(xf) BCflightangle4constraintPeAApA(xf, body, PeA, ApA, incT, AttR);
bcfun = @(xf) BCkeplerian3constraintPeAApA(xf, body, PeA, ApA, incT);
%bcfun = @(xf) BCkeplerian5constraintPeAApA(xf, body, PeA, ApA, incT, LANT, ArgPT);

% launchsite latitude
lat = deg2rad(28.608);
lng = 0;

[ rf vf xf indexes r_scale v_scale t_scale ] = launch(stages, body, bcfun, lat, lng, incT);

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
orbit.rf = norm(rf);
orbit.vf = norm(vf);
orbit.gamma = acosd( min(max(norm(cross(rf,vf))/(norm(rf)*norm(vf)), -1), 1) );
orbit.inclination = rad2deg(i);
orbit.LAN = rad2deg(O);
orbit.argument_of_periapsis = rad2deg(o);
orbit.true_anomaly = rad2deg(nu);
orbit.semi_latus_rectum = p / 1000;

disp(orbit);

