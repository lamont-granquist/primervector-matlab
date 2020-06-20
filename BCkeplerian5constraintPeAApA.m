function z = BCkeplerian5constraintPeAApA(xf, body, PeA, ApA, incT, LANT, ArgPT)
  PeR = PeA + body.r;
  ApR = ApA + body.r;
  smaT = ( PeR + ApR ) / 2;
  eccT = ( ApR - PeR ) / ( ApR + PeR );

  z = BCkeplerian5constraint(xf, smaT, eccT, incT, LANT, ArgPT);
end
