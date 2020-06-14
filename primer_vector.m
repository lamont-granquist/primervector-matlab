close all; clear all; clc;
format longG;

global g_bar r0_bar v0_bar r_scale v_scale

%
% CONSTANTS
%

mu     = 3.986004418e+14;  % m^3/s^2; earth
rearth = 6.371e+6;                % m; earth radius
wearth = 7.2921159e-5 * [0,0,1]'; % rad/s; sidereal angular

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
PeR = PeA + rearth;
ApR = ApA + rearth;
smaT = ( PeR + ApR ) / 2;
eccT = ( ApR - PeR ) / ( ApR + PeR );
rT = PeR;
vT = sqrt(mu * ( 2 / PeR - 1 / smaT ) );
gammaT = deg2rad(0);

bcfun = @(xf) BCflightangle4constraint(xf, incT, rT, vT, gammaT);

% launchsite latitude
lat = deg2rad(28.608);
lng = 0;

% setup launch site along x-axis rotated up to latitude
r0 = rearth * [ cos(lat)*cos(lng) cos(lat)*sin(lng) sin(lat) ]';
v0 = cross(wearth, r0);
dr = r0 / norm(r0);
de = cross([ 0, 0, 1 ]', dr);
de = de / norm(de);
dn = cross(dr, de);
dn = dn / norm(dn);

% guess the initial costate (heading from spherical trig, pitch of 45 degrees)
heading = asin( min(1, max(0, cos(incT)/ cos(lat))) );
dh = dn * cos(heading) + de * sin(heading);

pitch_guess = 45;

pv0 = dh * cosd(pitch_guess) + dr * sind(pitch_guess);
pr0 = r0 / norm(r0) * 8 / 3;

[rf vf] = primer_vector2(stages, bcfun, mu, r0, v0, pv0, pr0);

%
% SOME OUTPUT
%

[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(rf', vf', mu);

orbit.semi_major_axis = a / 1000;
orbit.eccentricity = eMag;
orbit.PeR = (1 - eMag) * a;
orbit.ApR = (1 + eMag) * a;
orbit.PeA = orbit.PeR - rearth;
orbit.ApA = orbit.ApR - rearth;
orbit.inclination = rad2deg(i);
orbit.LAN = rad2deg(O);
orbit.argument_of_periapsis = rad2deg(o);
orbit.true_anomaly = rad2deg(nu);
orbit.semi_latus_rectum = p / 1000;

disp(orbit);

function [ rf vf xf indexes r_scale v_scale ] = primer_vector2(stages, bcfun, mu, r0, v0, pv0, pr0)
  global indexes r_scale v_scale r0_bar v0_bar g_bar Nphases

  %
  % INTERNAL CONSTANTS
  %

  g0     = 9.80665;                 % m/s; standard gravity
  Nphases = length(stages);
  indexes.r = 1:3;
  indexes.v = 4:6;
  indexes.pv = 7:9;
  indexes.pr = 10:12;
  indexes.m = 13;
  indexes.integrated = 1:13;
  indexes.bt = 14;
  indexes.total = 14;

  %
  % PROBLEM SCALING
  %

  g_bar = mu / norm(r0)^2;
  r_scale = norm(r0);
  v_scale = sqrt( norm(r0) * g_bar );
  t_scale = sqrt( norm(r0) / g_bar );

  r0_bar = r0 / r_scale;
  v0_bar = v0 / v_scale;

  for p = 1:Nphases;
    stages(p).ve     = stages(p).isp * g0;
    stages(p).a0     = stages(p).thrust / stages(p).m0;
    stages(p).tau    = stages(p).ve / stages(p).a0;
    stages(p).c      = g0 * stages(p).isp / t_scale;
    stages(p).bt_bar = stages(p).bt / t_scale;
  end

  %
  % INITIAL STATE AND COSTATE GUESS
  %

  x0(indexes.r)  = r0_bar;
  x0(indexes.v)  = v0_bar;
  x0(indexes.pv) = pv0;
  x0(indexes.pr) = pr0;
  x0(indexes.m)  = stages(1).m0;

  %
  % MAIN SOLVER
  %

  x0 = singleShooting(x0, stages);
  options = optimset('Algorithm','levenberg-marquardt','TolX',1e-15,'TolFun',1e-15,'MaxFunEvals',20000,'MaxIter',3000);
  lb(1:length(x0)) = -inf;
  ub(1:length(x0)) = inf;
  [x, z, exitflag, output, jacobian] = lsqnonlin(@(x) residualFunction(x, stages, bcfun), x0, lb, ub, options);
  xf = multipleShooting(x, stages);

  rf = xf(indexes.r + indexes.total) * r_scale;
  vf = xf(indexes.v + indexes.total) * v_scale;
end

%
% RESIDUAL FUNCTION
%

function z = residualFunction(x0, stages, bcfun)
  global indexes r_scale v_scale r0_bar v0_bar g_bar Nphases

  % pin a few values the optimizer should not play with
  % (removing them completely would make for faster, but messier code)
  x0(indexes.r) = r0_bar;
  x0(indexes.v) = v0_bar;
  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    x0(indexes.m + i_offset) = stages(p).m0;  % FIXME: mass continuity (upper stage coast)
    x0(indexes.bt + i_offset) = stages(p).bt_bar;  % FIXME: variable burntime
  end

  xf = multipleShooting(x0, stages);

  z = [];

  % initial conditions
  z = [
    x0(indexes.r)' - r0_bar
    x0(indexes.v)' - v0_bar
    x0(indexes.m) - stages(1).m0
    ];

  % transversality conditions
  z = [
    z
    bcfun(xf(indexes.integrated + (Nphases-1)*indexes.total));
    ];

  % burntime of the bottom stage
  z = [
    z
    x0(indexes.bt)' - stages(1).bt_bar
    ];

  % initial conditions and burntime of the pth+1 phase (continuity, mass jettison and upper stage burntime)
  for p = 1:Nphases-1
    i_offset = (p-1)*indexes.total;
    i_offset2 = p*indexes.total;

    % continuity
    z = [
      z
      xf(indexes.r + i_offset)' - x0(indexes.r + i_offset2)'
      xf(indexes.v + i_offset)' - x0(indexes.v + i_offset2)'
      xf(indexes.pr + i_offset)' - x0(indexes.pr + i_offset2)'
      xf(indexes.pv + i_offset)' - x0(indexes.pv + i_offset2)'
      ];
    % mass jettison
    z = [
      z
      x0(indexes.m + i_offset2)' - stages(p+1).m0
      ];
    % burntime
    z = [
      z
      x0(indexes.bt + i_offset2)' - stages(p+1).bt_bar
      ];
  end
end

%
% SINGLE SHOOTING INITIALIZATION
%

function x0 = singleShooting(x0, stages)
  global indexes Nphases

  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    i_offset2 = p*indexes.total;
    bt = stages(p).bt_bar;
    x0(indexes.bt + i_offset) = bt;
    x0(indexes.m + i_offset) = stages(p).m0;
    index_range = indexes.integrated + i_offset;
    ode45options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [ts, xs] = ode45(@(t,x) EOM(t, x, p, stages), [0 bt], x0(index_range), ode45options);
    if p < Nphases
      x0(indexes.integrated + i_offset2) = xs(end,:);
    end
  end
end

%
% MULTIPLE SHOOTING INTEGRATION
%

function xf = multipleShooting(x0, stages)
  global indexes Nphases

  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    bt = x0(indexes.bt + i_offset);
    index_range = indexes.integrated + i_offset;
    ode45options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [ts, xs] = ode45(@(t,x) EOM(t, x, p, stages), [0 bt], x0(index_range), ode45options);
    xf(index_range) = xs(end,:);
  end
end

%
% EQUATIONS OF MOTION
%

function dX_dt = EOM(t, X, p, stages)
  global indexes g_bar

  thrust = stages(p).thrust;
  c      = stages(p).c;

  r  = X(indexes.r);
  v  = X(indexes.v);
  pv = X(indexes.pv);
  pr = X(indexes.pr);
  m  = X(indexes.m);

  u = pv/norm(pv);
  T = thrust / (m * g_bar);

  r2 = dot(r,r);
  r3 = r2^(3/2);
  r5 = r2 * r3;

  rdot  = v;
  vdot  = - r / r3 + T * u;
  pvdot = - pr;
  prdot = pv / r3 - 3 / r5 * dot(r, pv) * r;
  mdot  = - thrust / c;

  dX_dt = [ rdot' vdot' pvdot' prdot' mdot ]';
end
