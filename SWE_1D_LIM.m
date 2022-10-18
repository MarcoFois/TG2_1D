%% Shallow water equations for dam break problem
%{
   Equations: one dimensional shallow water 
                   
              dq/dt + div(F(q)) = S(q),    q = [h,hu], F = [hu, u^2 + .5gh^2]

   Numerical Method: Explicit two-step Second-order Taylor-Galerkin

   Flux Limiter: Lax-Friedrichs Flux 
   
%}

clear all;


%% Domain and mass matrix (lumped)

g  = 9.81;
x   = linspace (0, 20, 100) .';
dx  = diff (x);

mass = [dx(1)/2;dx(2:end);dx(end)/2];


T   = 1.5;
dt  = .005;
tsave = linspace (0, 100, 1001);

%% Initial conditions

 h0 = (x<8)*2 + (x>=8 )*1 ; 
u0 = zeros (size (x));

h = h0;
u = u0;

hsave(:,1) = h;
usave = u; 
epsilon = 5e-4;
ii = 1;
figure()
time = 0;
while (time < T)
    
    %% First step: low order solution computation
   
    Su = sourceu(h,u);
    F_h = fluxh(h,u); F_u = fluxu(h,u); 
    h_cell = (h(2:end)+h(1:end-1))/2.; 
    u_cell = (u(2:end)+u(1:end-1))/2.;
    Su_cell = (Su(2:end)+Su(1:end-1))/2.;
    
    spatial_der_h = (h(2:end)-h(1:end-1))./dx; 
    temp_der_h = (F_h(2:end)-F_h(1:end-1))./dx;
    
    spatial_der_u = (u(2:end)-u(1:end-1))./dx; 
    temp_der_u = (F_u(2:end)-F_u(1:end-1))./dx;
    
    h_rec = h_cell + spatial_der_h .* [-dx/2,dx/2] + temp_der_h .* dt;
    u_rec = u_cell + spatial_der_u .* [-dx/2,dx/2] + temp_der_u .* dt;
%     phi_cell_h = flux_limiter(h_rec, h_cell); phi_cell_u = flux_limiter(u_rec, u_cell);
%     phi_cell_h(:) = 0; phi_cell_u(:) = 0; 
    
    h_12 = h_cell - dt/2 .* temp_der_h;
    u_12 = u_cell - dt/2 .* temp_der_u + dt/2.*Su_cell;
    
    
     %% Second step: low order solution
    
    
    diff_coeff = sqrt(g*h) + abs(u./(h+epsilon));
    diff_coeff = .5*(diff_coeff(2:end)+diff_coeff(1:end-1));
    
   
    dF_h = .5 * (h(2:end)-h(1:end-1)).*diff_coeff; dF_u = .5 * (u(2:end)-u(1:end-1)).*diff_coeff;
    F_hl = [fluxh(h_cell(1),u_cell(1));fluxh(h_cell,u_cell)-dF_h;fluxh(h_cell(end),u_cell(end))]; 
    F_ul = [fluxu(h_cell(1),u_cell(1)) ;fluxu(h_cell,u_cell)-dF_u  ;fluxu(h_cell(end),u_cell(end))];
    
    h = h - dt./mass .* (F_hl(2:end)-F_hl(1:end-1));
    u = u - dt./mass .* (F_ul(2:end)-F_ul(1:end-1));
    
    %% Boundary conditions: Dirchlet
    
    u(1)= 0;
    u(end) = 0;
    
    
    %% Limiter: Lax-Friedrichs flux
    
    Snod = sourceu(h,u);
    
    F_ul = [fluxu(h_cell(1),u_cell(1));fluxu(h_cell,u_cell)-dF_u ;fluxu(h_cell(end),u_cell(end))];
    F_h = [fluxh(h_12(1),u_12(1));fluxh(h_12,u_12);fluxh(h_12(end),u_12(end))]-F_hl;
    F_u = [fluxu(h_12(1),u_12(1));fluxu(h_12,u_12);fluxu(h_12(end),u_12(end))]-F_ul;
    
    h_lim = flux_limiter(h, F_h);
    u_lim = flux_limiter(u, F_u);
    
    h = h + dt./mass .* h_lim.*(-F_h(2:end) + F_h(1:end-1));
    u = u + dt./mass .* u_lim.*(-F_u(2:end) + F_u(1:end-1)) + dt.*Snod;
    
    u(1)=0;
    u(end)=0;
    
    time = time + dt;

 phi = 0;
 
 %% Use this function to plot the dam-break exact solution
 [h_ex, Ux_ex] = exact_sol(phi, 8, time, x, 2, 1,0, 0);

 plot(x,h_ex,'-',x,h,'-x')
 axis([0 20 0.0 3])
 legend('exact','computed');
 title (sprintf ("t = %g", time))
 grid on;
 drawnow

hsave(:,ii) = h;
ii = ii+1;

end

return

%% Useful functions

function Fh = fluxh (h, u)
Fh = u;
end

function Sh = sourceh (h, u)
Sh = 0 * h;
end

function Su = sourceu (h, u)
tau = .009;
G = 9.81*sin(atan(0));
Su = G*h*0 ;
end

function Fu = fluxu (h, u)
g = 9.81;
epsilon = 1e-5;
gz = g*cos(atan(2))*0;
Fu = u.^2 ./ (h) + (1/2) * g* h.^2;
end

function lim = flux_limiter (h, F)
%toll=1e-7; 
hmin = min(min(h, [h(2:end);h(end)]), [h(1);h(1:end-1)]);
hmax = max(max(h, [h(2:end);h(end)]), [h(1);h(1:end-1)]);

F_mat = [-F(2:end), F(1:end-1)];
P_p = sum(max(0,F_mat),2);
P_m = sum(min(0,F_mat),2);

Q_p = (hmax-h);  
Q_m = (hmin-h); 

R_p = (P_p==0).*1 + (P_p~=0).*min(1,Q_p./(P_p+eps));
R_m = (P_m==0).*1 + (P_m~=0).*min(1,Q_m./(P_m+eps));
 
lim = min(R_p.*(F_mat>=0) + R_m.*(F_mat<0),[],2); 
  
end



