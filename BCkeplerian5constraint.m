function z = BCkeplerian5constraint(xf, smaT, eccT, incT, LANT, ArgPT)
  global indexes r_scale v_scale

  smaT_bar = smaT / r_scale;

  s1 = sin(LANT); s2 = sin(incT); s3 = sin(ArgPT);
  c1 = cos(LANT); c2 = cos(incT); c3 = cos(ArgPT);

  hT = [ s1*s2 -c1*s2 c2 ]' * sqrt( smaT_bar * (1-eccT*eccT) );
  eT = [ c3*c1-c2*s3*s1 c3*s1+c2*c1*s3 s3*s2 ]' * eccT;

  rf = xf(indexes.r)';
  vf = xf(indexes.v)';
  pvf = xf(indexes.pv)';
  prf = xf(indexes.pr)';

  if (abs(hT(3)) <= 1e-6)
    hT = circshift(hT, 1);
    eT = circshift(eT, 1);
    rf = circshift(rf, 1);
    vf = circshift(vf, 1);
    pvf = circshift(pvf, 1);
    prf = circshift(prf, 1);
  end

  hf = cross(rf, vf);
  ef = cross(vf, hf) - rf / norm(rf);

  rf3 = norm(rf) * dot(rf,rf);

  hmiss = hf - hT;
  emiss = ef - eT;

  z = [
    hmiss
    emiss(1:2)
    % 1 transversality condition
    dot(prf, vf) - dot(pvf, rf) / rf3
  ];
end
