% From lu 2012 "A Versatile Powered Guidance Algorithm"
% radius, velocity, flight path angle, inclination (PEG style constriants)
function z = BCflightangle4constraint(xf, incT, rT, vT, gammaT)
  global indexes r_scale v_scale

  rT_bar = rT / r_scale;
  vT_bar = vT / v_scale;

  rf = xf(indexes.r);
  vf = xf(indexes.v);
  pvf = xf(indexes.pv);
  prf = xf(indexes.pr);

  n_hat = [ 0, 0, 1 ];
  rn = cross(rf, n_hat);
  vn = cross(vf, n_hat);

  hf = cross(rf, vf);

  z = [
    % 4 constraints
    ( norm(rf)^2 - rT_bar^2 ) / 2
    ( norm(vf)^2 - vT_bar^2 ) / 2
    dot(n_hat, hf) - norm(hf) * cos(incT)
    dot(rf, vf) - norm(rf) * norm(vf) * sin(gammaT)
    % 2 transversality conditions
    rT_bar^2 * ( dot(vf, prf) - vT_bar * sin(gammaT) / rT_bar * dot(rf, prf) ) - vT_bar^2 * ( dot(rf, pvf) - rT_bar * sin(gammaT) / vT_bar * dot(vf, pvf) )
    dot(hf,prf) * dot(hf, rn) + dot(hf, pvf) * dot(hf, vn)
  ];
end
