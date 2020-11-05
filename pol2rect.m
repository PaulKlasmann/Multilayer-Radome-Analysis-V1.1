function [x, y] = pol2rect (ANG, MAG)
  
  x = MAG * cos(ANG);
  y = MAG * sin(ANG);
  
endfunction
