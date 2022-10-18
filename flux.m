function [Fh, FU] = flux (h, U, g)

  Fh = zeros(size(h)); 
  FU = zeros (size (h));
  v = U ./ h;
  v(h <= 0) = 0;
  
  Fh = U;
  FU = U .* v + g * (h .^2) / 2;   

end