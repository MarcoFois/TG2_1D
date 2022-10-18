%%  Shallow Waters - TG2-2S - Sloshing problem
clear all
clc
%% Domain and FE basis and gradient

L              = 400;
N              = 100;
msh.ndof       = N+1;
msh.nel        = msh.ndof - 1;
msh.x          = linspace (-L, L, msh.ndof).';
msh.xm         = (msh.x(2:end)+msh.x(1:end-1))/2;
msh.conn       = [1:msh.ndof-1; 2:msh.ndof];
msh.h          = diff (msh.x);
msh.shg(1, :)  = -1./msh.h;
msh.shg(2, :)  = +1./msh.h;
nit = 5;
damp = 1;

%% Mass
mass           =  sparse (diag (([msh.h;0] + [0;msh.h])/3) +...
                          diag (msh.h/6, 1) + diag (msh.h/6, -1));
Dmass          = diag (sum (mass, 2));
%%Dmass = diag (diag (mass));

%% Bed
gradz(1:msh.ndof,1) = 0;
%gradz(1:msh.ndof,1) = -2/(2*L);

g = 9.81;
Z = cumtrapz (msh.x, gradz);
Zm = (Z(2:end)+Z(1:end-1))/2;

%% Initial condition
h0 = 8 - sin (pi*msh.x/2/L) - Z;% 5 - Z;
U0 = zeros (size (msh.x));

h = h0;
U  = U0;

T    = 1500;
t    = 0;
it   = 1;
dt   = 1e-2;
Njac = 100;

nsave = round (T/5);
tsave = linspace (0, T, nsave);

hsave = h0;
Usave = U0;
  
for is = 2 : nsave 
  while (t(it) < tsave (is))

    t(it+1) = t(it) + dt;
    if (t(it+1) > tsave (is))
      t(it+1) = tsave (is);
      dt = t(it+1) - t(it);
    end
    fprintf ("it = %d, t = %g, dt = %g\n", it, t(it+1), dt)

   
    [Fh, FU] = flux (h, U, g);    
    hm = (h(2:end)+h(1:end-1))/2;
    Um = (U(2:end)+U(1:end-1))/2;
    [~,  SU] = source (h, U, gradz, g);

    dh = zeros (msh.nel, 1);
    dU = zeros (msh.nel, 1);
    for k = 1 : msh.nel
      dh(k) = dh(k) + ...
      -(Fh(msh.conn(2, k)) - Fh(msh.conn(1, k)));
    
      dU(k) = dU(k) + ...
      -(FU(msh.conn(2, k)) - FU(msh.conn(1, k))) + ...
      (SU(msh.conn(2, k)) + SU(msh.conn(1, k))) * msh.h(k)/2;
      
    end 
    
    %% First step
    wh = hm + (dt/2) * (dh./msh.h(:));
    wU = Um + (dt/2) * (dU./msh.h(:));

    dh = zeros (msh.ndof, 1);
    dU = zeros (msh.ndof, 1);
    [Fh, FU] = flux (wh, wU, g);
    [~,  SU] = source (wh, wU, (gradz(1:end-1)+gradz(2:end))/2, g);

    for k = 1 : msh.nel
      for i = 1 : 2
        dh(msh.conn(i, k)) = dh(msh.conn(i, k)) + ...
                             Fh(k) * msh.shg(i, k) * msh.h(k);

         dU(msh.conn(i, k)) = dU(msh.conn(i, k)) + ...
                              FU(k) * msh.shg(i, k) * msh.h(k) + ...
                              SU(k) * msh.h(k)/2;
      end
    end
    
    
    %% Second step
    h = h + dt * (jacobi (mass, Dmass, dh, Njac,damp));
    U(2:end-1) = U(2:end-1) + dt * (jacobi (mass(2:end-1,2:end-1), Dmass(2:end-1,2:end-1), dU(2:end-1), Njac,damp));
    
    
    U([1, end]) = 0;
    v = U ./ h;
    v(h <= 0) = 0;
    
    dtold = dt;
    dt = (.9/sqrt(3)) * min (min (msh.h) ./ (norm (v, inf) + norm (sqrt (g * h), inf)));
    dt = min (dt, 2 * dtold);        
    

it = it + 1;
    end
 
%   figure (3)
%   plot (msh.x, h+Z,msh.xm, wh+Zm , msh.x, Z , 'linewidth', .5)
%   grid on
%   axis ([-L, L, 6, 10])
%   title (sprintf ("t = %g", tsave(is)))
%   drawnow
  
    figure (3)
  plot (msh.x, h+Z,'b-',msh.x,Z,'linewidth',1.2)
  grid on
  axis ([-L, L, 0, 10])
  title (sprintf ("Time (s) = %g", tsave(is)))
  xlabel('x (m)');
  ylabel('height (m)');
  drawnow
  
  
  
%   figure (4)
%   plot (msh.x, U, 'b-', msh.xm, wU ,'r--', 'linewidth', 1.5)
%   title (sprintf ("t = %g", tsave(is)))
%   drawnow
  

  hsave (:, is) = h;
  Usave (:, is) = U;

  end