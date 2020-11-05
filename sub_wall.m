function [Tw, IPD, Rw] = sub_wall (theta, f, er, tand, Ko, t, N, pol)
  CTH = cos(theta);
  STH = sin(theta) * sin(theta);

  % Calculate complex permitivity for each layer
  for i = 1:N
    x = er(i);
    y = er(i) * -(tand(i));
    [ANG, MAG] = rect2pol (x, y);
    emag(i) = MAG;
    eph(i) = ANG;
  end
  % Calculate layer impedance for perpendicular polarisation
  for i = 1:N
    x = er(i) - STH;
    y = er(i) * (-tand(i));
    [ANG, MAG] = rect2pol (x, y);
    [RMAG, RANG] = complex_sqr_root(MAG, ANG);
    magterm(i) = RMAG;
    angterm(i) = RANG;
    phimag(i) = Ko * magterm(i);
    phiang(i) = angterm(i);
    MAG = phimag(i);
    ANG = phiang(i);
    [x, y] = pol2rect (ANG, MAG);
    rephi(i) = x;
    imphi(i) = y;
    zmag(i) = CTH / (magterm(i) + 1e-12);
    zph(i) = -angterm(i);
    MAG = zmag(i);
    ANG = zph(i);
    [x, y] = pol2rect (ANG, MAG);
    rez(i) = x;
    imz(i) = y;
    if (pol == 1)         % Calculate layer impedance for parallel polarisation
      zmag(i) = 1 / ((emag(i) * zmag(i)) + 1e-12);
      zph(i) = -(eph(i) + zph(i));
      MAG = zmag(i);
      ANG = zph(i);
      [x, y] = pol2rect (ANG, MAG);
      rez(i) = x;
      imz(i) = y;
    end  
  end
% End of layer impedance calculation
  zmag(N+1) = 1;
  zph(N+1) = 0;
  rez(N+1) = 1;
  imz(N+1) = 0;
% Do next loop i = 1 serperately from the loop to include z for space before radome
% because indexing starts at 1 and not 0. Then proceed with loop starting i=2
  i = 1;
  x = rez(i) - 1;
  y = imz(i) - 0;
  [ANG, MAG] = rect2pol (x, y);
  NUMMAG = MAG;
  NUMANG = ANG;
  x = rez(i) + 1;
  y = imz(i) + 0;
  [ANG, MAG] = rect2pol (x, y);
  DENMAG = MAG;
  DENANG = ANG;
  RMAG(i) = NUMMAG / DENMAG;
  Rph(i) = NUMANG - DENANG;
  MAG = RMAG(i);
  ANG = Rph(i);
  [x, y] = pol2rect (ANG, MAG);
  reR(i) = x;
  imR(i) = y;
  reT(i) = 1 + reR(i);
  imT(i) = imR(i);
  x = reT(i);
  y = imT(i);
  [ANG, MAG] = rect2pol (x, y);
  Tmag(i) = MAG;
  Tph(i) = ANG;

% Continue from i=2
  for i = 2:(N+1)
    x = rez(i) - rez(i-1);
    y = imz(i) - imz(i-1);
    [ANG, MAG] = rect2pol (x, y);
    NUMMAG = MAG;
    NUMANG = ANG;
    x = rez(i) + rez(i-1);
    y = imz(i) + imz(i-1);
    [ANG, MAG] = rect2pol (x, y);
    DENMAG = MAG;
    DENANG = ANG;
    RMAG(i) = NUMMAG / DENMAG;
    Rph(i) = NUMANG - DENANG;
    MAG = RMAG(i);
    ANG = Rph(i);
    [x, y] = pol2rect (ANG, MAG);
    reR(i) = x;
    imR(i) = y;
    reT(i) = 1 + reR(i);
    imT(i) = imR(i);
    x = reT(i);
    y = imT(i);
    [ANG, MAG] = rect2pol (x, y);
    Tmag(i) = MAG;
    Tph(i) = ANG;  
  end

% Matrix multiplications begin here
Amag(1) = exp(-imphi(1) * t(1));
Amag(4) = 1 / Amag(1);
Amag(2) = RMAG(1) * Amag(4);
Amag(3) = RMAG(1) * Amag(1);
Aph(1) = rephi(1) * t(1);
Aph(2) = Rph(1) - Aph(1);
Aph(3) = Rph(1) + Aph(1);
Aph(4) = -Aph(1);

for K = 2:N
  Bmag(1) = exp(-imphi(K) * t(K));
  Bmag(4) = 1 / Bmag(1);
  Bmag(2) = RMAG(K) * Bmag(4);
  Bmag(3) = RMAG(K) * Bmag(1);
  Bph(1) = rephi(K) * t(K);
  Bph(2) = Rph(K) - Bph(1);
  Bph(3) = Rph(K) + Bph(1);
  Bph(4) = -Bph(1);
  [Amag, Aph] = cmult(Amag, Aph, Bmag, Bph);
end

Bmag(1) = 1;
Bmag(4) = 1;
Bmag(2) = RMAG(N + 1);
Bmag(3) = Bmag(2);
Bph(1) = 0;
Bph(2) = Rph(N + 1);
Bph(3) = Bph(2);
Bph(4) = 0;
[Amag, Aph] = cmult(Amag, Aph, Bmag, Bph);
tranmag = 1;
tranph = 0;

for j = 1:N +1
  tranmag = tranmag * Tmag(j);
  tranph = tranph + Tph(j);
end

tranmag = 1 / (tranmag + 1e-12);

while (tranph > (2 * pi))
    tranph = tranph - (2 * pi);
endwhile

tranph = -tranph;
Tw = 1 / (tranmag * Amag(1));
Ttotph = -(tranph + Aph(1));
Rw = Amag(3) / Amag(1);
SUM = 0;
for j = 1:N
  SUM = SUM + t(j);  
end

% Get IPD into correct quadrant
do
  if(Aph(1) < (CTH * Ko * SUM))
    Aph(1) = Aph(1) + (2 * pi);
  end  
until (Aph(1) >= (CTH * Ko * SUM))

IPD = Aph(1) - (CTH * Ko * SUM);
endfunction