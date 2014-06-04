clear all
close all



Ae = 0.4;
At = 0.2;
x = linspace(0,1,1000);


A = (25/9)*(Ae-At)*(x-2/5).^2+At;


As = At;
gam = 1.4;
for i = 1:length(x)
   F = @(m) (A(i)/As)-(1/m)*((2/(gam+1))*(1+((gam-1)/2)*m^2))^((gam+1)/(2*(gam-1)));
   if (x(i) <= 2/5)
      M(i) = fzero(F,[0.01,1]); 
      
      M(i)
      x(i)
      F(M(i))
      
   else
      M(i) = fzero(F,[1,100]);
   end
    
end

rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
P = (1+((gam-1)/2)*M.^2).^(-gam/(gam-1));
T = (1+((gam-1)/2)*M.^2).^-1;
plot(x,M,x,A,x,P)

% syms x
% % syms a0 c1 c2 c3 cx1 cx2 cx3 ax1 ax2 ax3 gamma
% a0 = 2/3;
% c1 = 1;
% cx1 = .15;
% ax1 = 0.075*pi;
% c2 = 70;
% cx2 = 7;
% ax2 = 0.15*pi;
% c3 = 1e-5;
% cx3 = 2e-4;
% ax3 = 0.1*pi;
% gam = 1.4;
% 
% 
% 
% a = a0*(1+0.5*cos(2*pi*x));
% rho = c1+cx1*sin(ax1*x);
% u = c2+cx2*cos(ax2*x);
% p = c3+cx3*sin(ax3*x);
% 
% E = (1/(gam-1))*(p/rho)+0.5*u^2;
% 
% diff(rho*u*a,x);
% simplify(ans)
% 
% diff((rho*u^2+p)*a,x)-p*diff(a,x);
% simplify(ans)
% 
% diff(rho*u*(E+p/rho)*a,x);
% simplify(ans)
% 
