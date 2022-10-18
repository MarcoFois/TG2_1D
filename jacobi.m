function x = jacobi (A, P, b, nit, damp)
  
  res = b; 
  x   = damp*(P\ res);
  for it = 2 : nit
    x = x + damp*(P \ res);
    res = (b - A * x);
  end
  end