function [x, rho,u,P ] = initializeeuler(obj )
%INITIALIZEEULER Summary of this function goes here
%   Detailed explanation goes here
 Ae = 0.4;
 At = 0.2;
P0 = obj.P0;
T0 = obj.T0;
Pb = obj.Pb;%0.97*P0;

x = linspace(0,1,1000);
N = obj.nCells;
h = obj.cellWidths;
% A= (25/9)*(Ae-At)*(x-2/5).^2+At;
A = zeros(size(x));
for i = 1:length(x)
A(i) = obj.getArea(x(i));
end
As = At;
gam = 1.4;

fprintf('subsonic constant')
for i = 1:length(x)

%     F = @(m) (A(i)/As)-(1/m)*((2/(gam+1))*(1+((gam-1)/2)*m^2))^((gam+1)/(2*(gam-1)));
%    if ((x(i) <= 2/5 ))
%        
%       M(i) = fzero(F,[0.01,1]); 
%    else
%       M(i) = fzero(F,[1,100]);
%    end

M(i) = sqrt((2/(gam-1))*((P0/Pb)^((gam-1)/gam)-1));
end



% 
% rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
P = P0*(1+((gam-1)/2)*M.^2).^(-gam/(gam-1));

T = T0*(1+((gam-1)/2)*M.^2).^-1;
figure
subplot(3,1,1)
plot(x,M,x,A,x,P)
legend('M','A','P')

% [ sqrt(gam*P(1)/rho(1)) sqrt(T(1))]


rho = (1+((gam-1)/2)*M.^2).^(-1/(gam-1));
% T
% rho = P/T

u = sqrt(gam*P./rho).*M;
% error('1')


E = (1/(gam-1))*(P./rho)+0.5*u.^2;
% subplot(3,1,2)
% plot(x,rho,x,rho.*u,x,rho.*E)
% legend('\rho','\rho u','\rho E','Interpreter','Latex')
% subplot(3,1,3)
% plot(x,rho,x,u,x,P)
% legend('\rho','u','P','Interpreter','Latex')
% 



%%%
if(obj.areatype == 1)
As = Ae/( (1/M(end))*((2/(gam+1))*(1+((gam-1)/2)*M(end)^2))^((gam+1)/(2*(gam-1))) );
for i = 1:length(x)
    F = @(m) (A(i)/As)-(1/m)*((2/(gam+1))*(1+((gam-1)/2)*m^2))^((gam+1)/(2*(gam-1)));
  
      M(i) = fzero(F,[0.01,1]); 
  
    rho(i) = (1+((gam-1)/2)*M(i)^2)^(-1/(gam-1));
    P(i) = (1+((gam-1)/2)*M(i)^2)^(-gam/(gam-1));
    u(i) = M(i)*sqrt(gam*P(i)/rho(i));
end
end
%%%

rho(1)
rho(end)
u(1)
u(end)
P(1)
P(end)
% error('1')


subplot(3,1,2)
plot(x,rho,x,rho.*u,x,rho.*E)
legend('\rho','\rho u','\rho E','Interpreter','Latex')
subplot(3,1,3)
plot(x,rho,x,u,x,P)
legend('\rho','u','P','Interpreter','Latex')

% error('1')

rsp = spapi(8,x,rho);
rspi = fnint(rsp);
usp = spapi(8,x,u);
uspi = fnint(usp);
Psp = spapi(8,x,P);
Pspi = fnint(Psp);

xx = obj.cellCentroids;
vav = NaN*ones(N+2,3);
for i = 2:N+1
    vav(i,1) = (1/h(i))*(fnval(rspi,xx(i)+h(i)/2)-fnval(rspi,xx(i)-h(i)/2)); 
    vav(i,2) = (1/h(i))*(fnval(uspi,xx(i)+h(i)/2)-fnval(uspi,xx(i)-h(i)/2)); 
    vav(i,3) = (1/h(i))*(fnval(Pspi,xx(i)+h(i)/2)-fnval(Pspi,xx(i)-h(i)/2)); 
end
obj.exactSolutionV = vav;

fnval(rspi,xx(end)+h(end)/2)
fnplt(rspi)
% error('2')


rsp = spapi(8,x,rho);
rspi = fnint(rsp);
rusp = spapi(8,x,rho.*u);
ruspi = fnint(rusp);
rEsp = spapi(8,x,rho.*E);
rEspi = fnint(rEsp);

xx = obj.cellCentroids;
uav = NaN*ones(N+2,3);
for i = 2:N+1
    uav(i,1) = (1/h(i))*(fnval(rspi,xx(i)+h(i)/2)-fnval(rspi,xx(i)-h(i)/2)); 
    uav(i,2) = (1/h(i))*(fnval(ruspi,xx(i)+h(i)/2)-fnval(ruspi,xx(i)-h(i)/2)); 
    uav(i,3) = (1/h(i))*(fnval(rEspi,xx(i)+h(i)/2)-fnval(rEspi,xx(i)-h(i)/2)); 
end

obj.exactSolutionU = uav;

end

