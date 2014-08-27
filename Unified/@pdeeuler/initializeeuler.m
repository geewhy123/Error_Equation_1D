function [x, rho,u,P ] = initializeeuler(obj )
%INITIALIZEEULER Summary of this function goes here
%   Detailed explanation goes here
 

% new MMS stuff
a0 = 0.3;
a1 = 0.1;
c1 = 0.925;
cx1 = 0.025;
ax1 = 2*pi;

c2 = 0.25;
cx2 =0.25;
ax2 = pi;

c3 = 0.9;
cx3 = 0.07;
ax3 = 2*pi;

gam = obj.gamma;
N = obj.nCells;
V = NaN*ones(N+2,3);
f = NaN*ones(N+2,3);
x = obj.cellCentroids;
h = obj.cellWidths;

for i = 2:N+1
    xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;

V(i,1) = (1/h(i))* (c1*(xr-xl)+ (cx1/ax1)*(sin(ax1*xr)-sin(ax1*xl)));
V(i,2) = (1/h(i))* (c2*(xr-xl)+ (cx2/ax2)*(-cos(ax2*xr)+cos(ax2*xl)));
V(i,3) = (1/h(i))* (c3*(xr-xl)+ (cx3/ax3)*(sin(ax3*xr)-sin(ax3*xl)));

ar = obj.getArea(xr);
al = obj.getArea(xl);
rhor = (c1+cx1*cos(ax1*xr));
rhol = (c1+cx1*cos(ax1*xl));
ur = (c2+cx2*sin(ax2*xr));
ul = (c2+cx2*sin(ax2*xl));
Pr = (c3+cx3*cos(ax3*xr));
Pl = (c3+cx3*cos(ax3*xl));
Er = (1/(gam-1))*Pr/rhor + 0.5*ur^2;
El = (1/(gam-1))*Pl/rhol + 0.5*ul^2;

f(i,1) = (1/h(i))*( rhor*ur*ar - rhol*ul*al ) ;
f(i,2) = (1/h(i))*( (rhor*ur^2+Pr)*ar - (rhol*ul^2+Pl)*al ) - (1/h(i))*( (a1*cos(2*pi*xr)*(2*c3 + cx3*cos(2*pi*xr)))/2 - (a1*cos(2*pi*xl)*(2*c3 + cx3*cos(2*pi*xl)))/2);
f(i,3) = (1/h(i))*( (rhor*ur*(Er+Pr/rhor))*ar -(rhol*ul*(El+Pl/rhol))*al );
end
assert(obj.areatype == 2);
%     f
%     g = pi*cx2*cos(pi*x).*(c1+cx1*cos(2*pi*x))-2*pi*cx1*sin(2*pi*x).*(c2+cx2*sin(pi*x))
%     u0 = c2;
%     u1 = cx2;
%     rho0 = c1;
%     rho1 = cx1;
%     P0 = c3;
%     P1 = cx3;
%     g = 2*pi*u1*cos(pi*x).*(rho0 + rho1*cos(2*pi*x)).*(u0 + u1*sin(pi*x)) - 2*pi*rho1*sin(2*pi*x).*(u0 + u1*sin(pi*x)).^2 - 2*P1*pi*sin(2*pi*x) - (P0+P1*cos(2*pi*x)).*(-2*pi*a1*sin(2*pi*x))
% %     g =(pi*(6*rho0*u1^3*cos(pi*x) - 4*rho1*u1^3*cos(pi*x) - 6*rho0*u1^3*cos(3*pi*x) + 9*rho1*u1^3*cos(3*pi*x) - 5*rho1*u1^3*cos(5*pi*x) - 16*rho1*u0^3*sin(2*pi*x) + 56*P0*u1*cos(pi*x) - 28*P1*u1*cos(pi*x) + 84*P1*u1*cos(3*pi*x) - 112*P1*u0*sin(2*pi*x) + 24*rho0*u0^2*u1*cos(pi*x) - 12*rho1*u0^2*u1*cos(pi*x) + 36*rho1*u0^2*u1*cos(3*pi*x) + 24*rho0*u0*u1^2*sin(2*pi*x) - 24*rho1*u0*u1^2*sin(2*pi*x) + 24*rho1*u0*u1^2*sin(4*pi*x)))/16
%     max(abs(f(:,2)-g))
%     figure
%     plot(x,f(:,2),x,g)
%     error('1')
f
% error('1')
obj.exactSolutionV = V;
obj.source = f;
x = linspace(0,1,1000);
rho = c1 + cx1*cos(ax1*x);
u = c2 +cx2*sin(ax2*x);
P = c3+cx3*cos(ax3*x);
% obj.exactSolutionU = U;
return;




% Ae = 0.4;
%  At = 0.2;
% P0 = obj.P0;
% T0 = obj.T0;
% Pb = obj.Pb;%0.97*P0;
% 
% x = linspace(0,1,1000);
% N = obj.nCells;
% h = obj.cellWidths;
% % A= (25/9)*(Ae-At)*(x-2/5).^2+At;
% A = zeros(size(x));
% for i = 1:length(x)
% A(i) = obj.getArea(x(i));
% end
% As = At;
% gam = 1.4;
% 
% fprintf('Quasi-1D Euler subsonic \n')
% for i = 1:length(x)
% 
% %     F = @(m) (A(i)/As)-(1/m)*((2/(gam+1))*(1+((gam-1)/2)*m^2))^((gam+1)/(2*(gam-1)));
% %    if ((x(i) <= 2/5 ))
% %        
% %       M(i) = fzero(F,[0.01,1]); 
% %    else
% %       M(i) = fzero(F,[1,100]);
% %    end
% 
% M(i) = sqrt((2/(gam-1))*((P0/Pb)^((gam-1)/gam)-1));
% end
% 
% 
% 
% % 
% % rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
% P = P0*(1+((gam-1)/2)*M.^2).^(-gam/(gam-1));
% 
% T = T0*(1+((gam-1)/2)*M.^2).^-1;
% % figure
% % subplot(3,1,1)
% % plot(x,M,x,A,x,P)
% % legend('M','A','P')
% 
% % [ sqrt(gam*P(1)/rho(1)) sqrt(T(1))]
% 
% 
% rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
% % T
% % rho = P/T
% 
% u = sqrt(gam*P./rho).*M;
% % error('1')
% 
% 
% E = (1/(gam-1))*(P./rho)+0.5*u.^2;
% % subplot(3,1,2)
% % plot(x,rho,x,rho.*u,x,rho.*E)
% % legend('\rho','\rho u','\rho E','Interpreter','Latex')
% % subplot(3,1,3)
% % plot(x,rho,x,u,x,P)
% % legend('\rho','u','P','Interpreter','Latex')
% % 
% 
% 
% 
% %%%
% if(obj.areatype == 1)
% As = Ae/( (1/M(end))*((2/(gam+1))*(1+((gam-1)/2)*M(end)^2))^((gam+1)/(2*(gam-1))) );
% for i = 1:length(x)
%     F = @(m) (A(i)/As)-(1/m)*((2/(gam+1))*(1+((gam-1)/2)*m^2))^((gam+1)/(2*(gam-1)));
%   
%       M(i) = fzero(F,[0.01,1]); 
%   
%     rho(i) = (1+((gam-1)/2)*M(i)^2)^(-1/(gam-1));
%     P(i) = (1+((gam-1)/2)*M(i)^2)^(-gam/(gam-1));
%     u(i) = M(i)*sqrt(gam*P(i)/rho(i));
% end
% end
% %%%
% 
% E = (1/(gam-1))*(P./rho)+0.5*u.^2;
% rho(1);
% rho(end);
% u(1);
% u(end);
% P(1);
% P(end);
% % error('1')
% figure
% subplot(3,1,1)
% plot(x,M,x,A,x,P)
% legend('M','A','P')
% 
% subplot(3,1,2)
% plot(x,rho,x,rho.*u,x,rho.*E);
% l1 = legend('$\rho$','$\rho u$','$\rho E$');
% set(l1,'Interpreter','Latex')
% subplot(3,1,3)
% plot(x,rho,x,u,x,P)
% l2 = legend('$\rho$','u','P');
% set(l2,'Interpreter','Latex')
% 
% % error('1')
% 
% rsp = spapi(8,x,rho);
% rspi = fnint(rsp);
% usp = spapi(8,x,u);
% uspi = fnint(usp);
% Psp = spapi(8,x,P);
% Pspi = fnint(Psp);
% 
% xx = obj.cellCentroids;
% vav = NaN*ones(N+2,3);
% for i = 2:N+1
%     vav(i,1) = (1/h(i))*(fnval(rspi,xx(i)+h(i)/2)-fnval(rspi,xx(i)-h(i)/2)); 
%     vav(i,2) = (1/h(i))*(fnval(uspi,xx(i)+h(i)/2)-fnval(uspi,xx(i)-h(i)/2)); 
%     vav(i,3) = (1/h(i))*(fnval(Pspi,xx(i)+h(i)/2)-fnval(Pspi,xx(i)-h(i)/2)); 
% end
% obj.exactSolutionV = vav;
% 
% fnval(rspi,xx(end)+h(end)/2);
% figure
% subplot(3,1,1)
% fnplt(rsp)
% % hold on
% % plot(xx,c1+cx1*cos(ax1*xx))
% subplot(3,1,2)
% fnplt(usp)
% % hold on
% % plot(xx,c2+cx2*sin(ax2*xx))
% subplot(3,1,3)
% fnplt(Psp)
% % hold on
% % plot(xx,c3+cx3*cos(ax3*xx))
% % error('2')
% 
% 
% rsp = spapi(8,x,rho);
% rspi = fnint(rsp);
% rusp = spapi(8,x,rho.*u);
% ruspi = fnint(rusp);
% rEsp = spapi(8,x,rho.*E);
% rEspi = fnint(rEsp);
% 
% xx = obj.cellCentroids;
% uav = NaN*ones(N+2,3);
% for i = 2:N+1
%     uav(i,1) = (1/h(i))*(fnval(rspi,xx(i)+h(i)/2)-fnval(rspi,xx(i)-h(i)/2)); 
%     uav(i,2) = (1/h(i))*(fnval(ruspi,xx(i)+h(i)/2)-fnval(ruspi,xx(i)-h(i)/2)); 
%     uav(i,3) = (1/h(i))*(fnval(rEspi,xx(i)+h(i)/2)-fnval(rEspi,xx(i)-h(i)/2)); 
% end
% 
% obj.exactSolutionU = uav;

end

