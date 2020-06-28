function [ rf vf xf indexes r_scale v_scale t_scale ] = primer_vector(phases, body, bcfun, r0, v0, pv0, pr0)
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
  indexes.p = 7:12;
  indexes.m = 13;
  indexes.integrated = 1:13;
  indexes.bt = 14;
  indexes.total = 14;

  %
  % PROBLEM SCALING / SANITY
  %

  g_bar = body.mu / norm(r0)^2;
  r_scale = norm(r0);
  v_scale = sqrt( norm(r0) * g_bar );
  t_scale = sqrt( norm(r0) / g_bar );

  r0_bar = r0 / r_scale;
  v0_bar = v0 / v_scale;

  for p = 1:Nphases;
    if ~isfield(phases(p), 'infinte') || isempty(phases(p).infinite)
      phases(p).infinite = false;
    end
    if ~isfield(phases(p), 'bt_free') || isempty(phases(p).bt_free)
      phases(p).bt_free = false;
    end
    if ~isfield(phases(p), 'unguided') || isempty(phases(p).unguided)
      phases(p).unguided = false;
    end
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
  options = optimset('Algorithm','levenberg-marquardt','TolX',1e-15,'TolFun',1e-15,'MaxFunEvals',2000,'MaxIter',300);
  lb(1:length(x0)) = -inf;
  ub(1:length(x0)) = inf;
  [x, z, exitflag, output, jacobian] = lsqnonlin(@(x) residualFunction(x, phases, bcfun), x0, lb, ub, options)
  xf = multipleShooting(x, phases);

  rf = xf(indexes.r + (Nphases-1) * indexes.total) * r_scale;
  vf = xf(indexes.v + (Nphases-1) * indexes.total) * v_scale;
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
    if ~phases(p).bt_free
      x0(indexes.bt + i_offset) = phases(p).bt_bar;
    end
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
  z = vertcat(z, burntimeResidual(phases, 1, x0, xf));

  % initial conditions and burntime of the pth+1 phase (continuity, mass jettison and upper stage burntime)
  for p = 2:Nphases
    i_offset = (p-2)*indexes.total;
    i_offset2 = (p-1)*indexes.total;

    % continuity
    z = [
      z
      xf(indexes.r + i_offset)' - x0(indexes.r + i_offset2)'
      xf(indexes.v + i_offset)' - x0(indexes.v + i_offset2)'
      xf(indexes.pr + i_offset)' - x0(indexes.pr + i_offset2)'
      xf(indexes.pv + i_offset)' - x0(indexes.pv + i_offset2)'
      ];
    % mass jettison or continuity
    z = [
      z
      x0(indexes.m + i_offset2)' - phases(p).m0
      ];
    % burntime of p+1th phase
    z = vertcat(z, burntimeResidual(phases, p, x0, xf));
  end
end

function z = burntimeResidual(phases, p, x0, xf)
  global Nphases indexes
  i_offset = (p-1)*indexes.total;

  r0  = x0(indexes.r + i_offset);
  v0  = x0(indexes.v + i_offset);
  pv0 = x0(indexes.pv + i_offset);
  pr0 = x0(indexes.pr + i_offset);
  H0t0 = dot(pr0, v0) - dot(pv0, r0) / norm(r0)^(3/2);

  rf  = xf(indexes.r + i_offset);
  vf  = xf(indexes.v + i_offset);
  pvf = xf(indexes.pv + i_offset);
  prf = xf(indexes.pr + i_offset);
  H0tf = dot(prf, vf) - dot(pvf, rf) / norm(rf)^(3/2);

  if phases(p).bt_free
    % free final time constraint on the mangnitude of the total costate vector for the optimized burntime
    % QUESTION: why does this seem to work even if we're optimizing the time of a prior stage?
    z = norm(xf(indexes.p + i_offset)) - 1;
  else
    % fixed burntime
    z = x0(indexes.bt + i_offset) - phases(p).bt_bar;
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
    phases(p).x0 = x0(index_range);
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
    phases(p).x0 = x0(index_range);
    [ts, xs] = ode45(@(t,x) EOM(t, x, p, phases), [0 bt], x0(index_range), ode45options);
    xf(index_range) = xs(end,:);
    xf(indexes.bt + i_offset) = bt;
  end
end

%
% EQUATIONS OF MOTION
%

function dX_dt = EOM(t, X, p, phases)
  global indexes g_bar

  thrust   = phases(p).thrust;
  c        = phases(p).c;
  infinite = phases(p).infinite;

  % X is only the integrated state and costate from a single phase, since those values come first
  % in any phase, followed by the time parameter, we can re-use the same indexes.
  r  = X(indexes.r);
  v  = X(indexes.v);
  pv = X(indexes.pv);
  pr = X(indexes.pr);
  m  = X(indexes.m);

  if phases(p).unguided
    % for inertially-fixed upper solid motors -- this constraint does not affect the costate equations, i think?
    pv0 = phases(p).x0(indexes.pv)';
    u = pv0 / norm(pv0);
  else
    u = pv/norm(pv);
  end
  T = thrust / (m * g_bar);

  if infinite; T = T * 2; end

  r2 = dot(r,r);
  r3 = r2^(3/2);
  r5 = r2 * r3;

  rdot  = v;
  vdot  = - r / r3 + T * u;
  pvdot = - pr;
  prdot = pv / r3 - 3 / r5 * dot(r, pv) * r;

  if ( infinite ) || ( thrust == 0 )
    mdot = 0;
  else
    mdot  = - thrust / c;
  end

  dX_dt = [ rdot' vdot' pvdot' prdot' mdot ]';
end
