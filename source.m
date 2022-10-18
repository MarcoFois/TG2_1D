function [Sh, SU] = source (h, U, gradz, g)

  Sh = zeros(size(h));
  SU = zeros (size (h));
  v = U ./ h;
  v(h <= 0) = 0;
  
  Sh(:)   = 0;
  SU      = -g*h.*gradz;
  
end