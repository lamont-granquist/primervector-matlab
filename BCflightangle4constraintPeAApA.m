function z = BCflightangle4constraintPeAApA(xf, body, PeA, ApA, incT, AttA);
  PeR  = PeA + body.r;
  ApR  = ApA + body.r;
  AttR = AttA + body.r;
  smaT = ( PeR + ApR ) / 2;
  eccT = ( ApR - PeR ) / ( ApR + PeR );
  hT = sqrt( body.mu * smaT * ( 1 - eccT * eccT ) );

  AttR = min(max(AttR, PeR), ApR);
  vT = sqrt( body.mu * ( 2/AttR - 1/smaT ) );
  gammaT = acos( min(max(hT / (AttR * vT), -1), 1) );

  z = BCflightangle4constraint(xf, AttR, vT, incT, gammaT);
end
