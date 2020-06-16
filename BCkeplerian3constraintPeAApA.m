function z = BCkeplerian3constraintPeAApA(xf, body, PeA, ApA, incT)
  PeR = PeA + body.r;
  ApR = ApA + body.r;
  smaT = ( PeR + ApR ) / 2;
  eccT = ( ApR - PeR ) / ( ApR + PeR );

  z = BCkeplerian3constraint(xf, smaT, eccT, incT);
end
