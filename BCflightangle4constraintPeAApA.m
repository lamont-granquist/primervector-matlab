function z = BCflightangle4constraintPeAApA(xf, body, PeA, ApA, incT, gammaT);
  PeR = PeA + body.r;
  ApR = ApA + body.r;
  smaT = ( PeR + ApR ) / 2;
  eccT = ( ApR - PeR ) / ( ApR + PeR );

  % FIXME: i'm just ignoring gammaT here
  rT = PeR;
  vT = sqrt(body.mu * ( 2 / PeR - 1 / smaT ) );

  z = BCflightangle4constraint(xf, rT, vT, incT, gammaT);
end
