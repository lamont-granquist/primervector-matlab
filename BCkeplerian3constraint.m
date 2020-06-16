% The transversality conditions for free LAN, ArgP and attachment are from:
%
% Beinfeng Pan, "Reduced Transversality Conditions in Optimal Space Trajectories"
% https://www.researchgate.net/publication/260671510_Reduced_Transversality_Conditions_in_Optimal_Space_Trajectories
%
function z = BCkeplerian3constraint(xf, smaT, eccT, incT)
  global indexes r_scale v_scale

  smaT_bar = smaT / r_scale;
  hTm = sqrt( smaT_bar * ( 1 - eccT * eccT ) );

  rf = xf(indexes.r);
  vf = xf(indexes.v);
  pvf = xf(indexes.pv);
  prf = xf(indexes.pr);

  n = [ 0, 0, 1 ];
  rn = cross(rf, n);
  vn = cross(vf, n);

  hf = cross(rf, vf);
  smaf = 1 / ( 2 / norm(rf) - dot(vf,vf));
  eccf = sqrt(1 - dot(hf,hf)/smaf);

  rf3 = norm(rf) * dot(rf,rf);

  z = [
    dot(hf,hf)/2 - hTm^2/2                              % angular momentum
    smaT_bar * ( 1 - eccT ) - ( smaf * ( 1 - eccf ) )   % periapsis constraint
    dot(n, hf / norm(hf)) - cos(incT)                   % inclination constraint
    %1 / ( 2 * smaf) - 1 / ( 2 * smaT_bar )            % energy constraint
    % eccT - eccf                                   % eccentricity constraint
    % eccT^2 - eccf^2                               % eceentricity^2 constraint
    % smaT_bar - smaf                                   % sma constraint
    % 3 transversality conditions
    dot(cross(prf, rf) + cross(pvf, vf), hf)
    dot(cross(prf, rf) + cross(pvf, vf), n)
    dot(prf, vf) - dot(pvf, rf) / rf3
  ];
end
