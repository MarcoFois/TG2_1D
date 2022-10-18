function[h_ex, Ux_ex] = exact_sol(phi, dam, t, X, hl, hr, ul, ur)
%% Function to plot the dam-break exact solution
%{

    INPUT:
        phi --------> bed slope (rad)
        dam --------> dam position on the domain
        t   --------> time (s)
        hl  --------> left height (m)
        hr  --------> right height (m)
        ul  --------> left mass flux (m^2/s^2)
        ur  --------> right mass flux (m^2/s^2)

    OUTPUT:
        h_ex -------> exact solution for h
        Ux_ex ------> exact solution for hu

%}


g = 9.81;
gz = g*cos(phi);
m  = g*sin(phi); % checked



nf = @(x) 2*(sqrt(gz*x)-sqrt(gz*hl)) + sqrt(.5*gz*(x+hr)/x/(hr+eps))*(x-hr)+ur-ul;
h_st = fzero(nf,(hr+hl)/2); % checked

u_st = 2*(sqrt(gz*hl)-sqrt(gz*h_st))+ul;


c_st = sqrt(gz*h_st);

cl = sqrt(gz*hl); cr = sqrt(gz*hr);
 s_dot_2 = ur + h_st*sqrt(gz/2*(h_st+hr)/(h_st*hr)); 




h_ex = X*0;
h_ex(X <= +dam + m*t^2/2+(ul-cl)*t) = hl;
h_ex(X > + dam + m*t^2/2+s_dot_2*t) = hr;
h_ex(X < dam + m*t^2/2+s_dot_2*t & X > dam +(ul+cl*2-3*c_st)*t+m*t^2/2) = h_st;

X_rar = X(X <= dam +(m*t^2/2+(ul+cl*2-3*c_st)*t) & X > dam +m*t^2/2+(ul-cl)*t);
h_ex(X <= (dam +m*t^2/2+(ul+cl*2-3*c_st)*t) & X >dam + m*t^2/2+(ul-cl)*t) = (ul+cl*2-((X_rar-dam)./t-m*t/2)).^2/9/gz;



Ux_ex = X*0;
Ux_ex(X <=dam + m*t^2/2+(ul-cl)*t) = hl*(ul+m*t);
Ux_ex(X > dam + m*t^2/2+s_dot_2*t) = hr*(ur+m*t);
Ux_ex(X < dam + m*t^2/2+s_dot_2*t & X >dam +(ul+cl*2-3*c_st)*t+m*t^2/2) = h_st*(u_st+m*t);

X_rar = X(X <= dam +(m*t^2/2+(ul+cl*2-3*c_st)*t) & X >dam + m*t^2/2+(ul-cl)*t);
Ux_ex(X <= (dam +m*t^2/2+(ul+cl*2-3*c_st)*t) & X >dam + m*t^2/2+(ul-cl)*t) = (ul+cl*2-((X_rar-dam)./t-m*t/2)).^2/9/gz .* ((ul+cl*2+2*((X_rar-dam)./t-m*t/2))/3+m*t);



end