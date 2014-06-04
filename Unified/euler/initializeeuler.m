function [ output_args ] = initializeeuler(obj )
%INITIALIZEEULER Summary of this function goes here
%   Detailed explanation goes here
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
%       
%       M(i)
%       x(i)
%       F(M(i))
   else
      M(i) = fzero(F,[1,100]);
   end
    
end

rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
P = (1+((gam-1)/2)*M.^2).^(-gam/(gam-1));
T = (1+((gam-1)/2)*M.^2).^-1;
figure
subplot(3,1,1)
plot(x,M,x,A,x,P)
legend('M','A','P')
u = sqrt(gam*P./rho).*M;
E = (1/(gam-1))*(P./rho)+0.5*u.^2;
subplot(3,1,2)
plot(x,rho,x,rho.*u,x,rho.*E)
legend('\rho','\rho u','\rho E','Interpreter','Latex')
subplot(3,1,3)
plot(x,rho,x,u,x,P)
legend('\rho','u','P','Interpreter','Latex')



rho(1)
rho(end)
u(1)
u(end)
P(1)
P(end)
end

