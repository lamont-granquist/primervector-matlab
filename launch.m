function [ rf vf xf indexes r_scale v_scale t_scale ] = launch(stages, body, bcfun, lat, lng, incT)
  % setup launch site along x-axis rotated up to latitude
  r0 = body.r * [ cos(lat)*cos(lng) cos(lat)*sin(lng) sin(lat) ]';
  v0 = cross(body.w, r0);
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

  phases = stages;

  stages(2)

  upper_burntime = -1;
  while upper_burntime < 0
    phases(length(phases)).infinite = true;

    [ rf vf xf indexes r_scale v_scale t_scale ] = primer_vector(phases, body, bcfun, r0, v0, pv0, pr0);

    upper_burntime = xf(indexes.bt + length(phases)-1*indexes.total);

    if upper_burntime < 0
      phases = phases(1:length(phases)-1);
    end
  end

  phases(length(phases)).infinite = false;

  [ rf vf xf indexes r_scale v_scale t_scale ] = primer_vector(phases, body, bcfun, r0, v0, pv0, pr0);

end

