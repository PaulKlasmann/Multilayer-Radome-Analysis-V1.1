function [Amag, Aph] = cmult(Amag, Aph, Bmag, Bph)
  for i = 1:2:3
    for j = 1:2
      R1 = Amag(i) * Bmag(j);
      R2 = Amag(i + 1) * Bmag(j + 2);
      E1 = Aph(i) + Bph(j);
      E2 = Aph(i + 1) + Bph(j + 2);
      R3 = (R1 * cos(E1)) + (R2 * cos(E2));
      E3 = (R1 * sin(E1)) + (R2 * sin(E2));
      Cmag(j + i - 1) = sqrt((R3*R3) + (E3*E3));
      if((R3 > 0) && (E3 > 0))
        PA = 0;
      elseif(R3 < 0)
        PA = pi;
      elseif((R3 > 0) && (E3 < 0))
        PA = 2 * pi;
      end
      Cph(j + i - 1) = atan(E3 / (R3 + 1e-12)) + PA;
    end
  end
  
  for i = 1:4
    Amag(i) = Cmag(i);
    Aph(i) = Cph(i);
  end
  
endfunction
