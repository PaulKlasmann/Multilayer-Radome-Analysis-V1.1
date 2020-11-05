function [ANG, MAG] = rect2pol (x, y)
  
  MAG = sqrt((x^2) + (y^2));
  
  if ((x>0) && (y>0))
    corr = 0;
  elseif ((x<0) && (y<0))
    corr = pi;
  elseif ((x>0) && (y<0))
    corr = 0;
  elseif (x<0) && (y>0)
    corr = pi;
  elseif (x>0) && (y==0)
    corr = 0;
  elseif (x<0) && (y==0)
    corr = pi;
  elseif (x==0) && (y!=0)
    corr = 0;
  elseif (x==0) && (y==0)
    ANG = 0;
    return;
  end
  
  ANG = atan(y / (x + 1e-12)) + corr;
  
endfunction
