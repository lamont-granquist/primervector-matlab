function [ rf vf xf indexes r_scale v_scale ] = primer_vector(phases, body, bcfun, r0, v0, pv0, pr0)
  global indexes r_scale v_scale r0_bar v0_bar g_bar Nphases

  %
  % INTERNAL CONSTANTS
  %

  g0     = 9.80665;                 % m/s; standard gravity
  Nphases = length(phases);
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

  g_bar = body.mu / norm(r0)^2;
  r_scale = norm(r0);
  v_scale = sqrt( norm(r0) * g_bar );
  t_scale = sqrt( norm(r0) / g_bar );

  r0_bar = r0 / r_scale;
  v0_bar = v0 / v_scale;

  for p = 1:Nphases;
    phases(p).ve     = phases(p).isp * g0;
    phases(p).a0     = phases(p).thrust / phases(p).m0;
    phases(p).tau    = phases(p).ve / phases(p).a0;
    phases(p).c      = g0 * phases(p).isp / t_scale;
    phases(p).bt_bar = phases(p).bt / t_scale;
  end

  %
  % INITIAL STATE AND COSTATE GUESS
  %

  x0(indexes.r)  = r0_bar;
  x0(indexes.v)  = v0_bar;
  x0(indexes.pv) = pv0;
  x0(indexes.pr) = pr0;
  x0(indexes.m)  = phases(1).m0;

  %
  % MAIN SOLVER
  %

  x0 = singleShooting(x0, phases);
  options = optimset('Algorithm','levenberg-marquardt','TolX',1e-15,'TolFun',1e-15,'MaxFunEvals',20000,'MaxIter',3000);
  lb(1:length(x0)) = -inf;
  ub(1:length(x0)) = inf;
  [x, z, exitflag, output, jacobian] = lsqnonlin(@(x) residualFunction(x, phases, bcfun), x0, lb, ub, options);
  xf = multipleShooting(x, phases);

  rf = xf(indexes.r + indexes.total) * r_scale;
  vf = xf(indexes.v + indexes.total) * v_scale;
end

%
% RESIDUAL FUNCTION
%

function z = residualFunction(x0, phases, bcfun)
  global indexes r_scale v_scale r0_bar v0_bar g_bar Nphases

  % pin a few values the optimizer should not play with
  % (removing them completely would make for faster, but messier code)
  x0(indexes.r) = r0_bar;
  x0(indexes.v) = v0_bar;
  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    x0(indexes.m + i_offset) = phases(p).m0;  % FIXME: mass continuity (upper stage coast)
    x0(indexes.bt + i_offset) = phases(p).bt_bar;  % FIXME: variable burntime
  end

  xf = multipleShooting(x0, phases);

  z = [];

  % initial conditions
  z = [
    x0(indexes.r)' - r0_bar
    x0(indexes.v)' - v0_bar
    x0(indexes.m) - phases(1).m0
    ];

  % transversality conditions
  z = [
    z
    bcfun(xf(indexes.integrated + (Nphases-1)*indexes.total));
    ];

  % burntime of the bottom stage
  z = [
    z
    x0(indexes.bt)' - phases(1).bt_bar
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
      x0(indexes.m + i_offset2)' - phases(p+1).m0
      ];
    % burntime
    z = [
      z
      x0(indexes.bt + i_offset2)' - phases(p+1).bt_bar
      ];
  end
end

%
% SINGLE SHOOTING INITIALIZATION
%

function x0 = singleShooting(x0, phases)
  global indexes Nphases

  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    i_offset2 = p*indexes.total;
    bt = phases(p).bt_bar;
    x0(indexes.bt + i_offset) = bt;
    x0(indexes.m + i_offset) = phases(p).m0;
    index_range = indexes.integrated + i_offset;
    ode45options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [ts, xs] = ode45(@(t,x) EOM(t, x, p, phases), [0 bt], x0(index_range), ode45options);
    if p < Nphases
      x0(indexes.integrated + i_offset2) = xs(end,:);
    end
  end
end

%
% MULTIPLE SHOOTING INTEGRATION
%

function xf = multipleShooting(x0, phases)
  global indexes Nphases

  for p = 1:Nphases
    i_offset = (p-1)*indexes.total;
    bt = x0(indexes.bt + i_offset);
    index_range = indexes.integrated + i_offset;
    ode45options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [ts, xs] = ode45(@(t,x) EOM(t, x, p, phases), [0 bt], x0(index_range), ode45options);
    xf(index_range) = xs(end,:);
  end
end

%
% EQUATIONS OF MOTION
%

function dX_dt = EOM(t, X, p, phases)
  global indexes g_bar

  thrust = phases(p).thrust;
  c      = phases(p).c;

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
