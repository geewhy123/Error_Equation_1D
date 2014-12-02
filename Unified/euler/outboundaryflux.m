function [ F1,F2,F3 ] = outboundaryflux( U1,U2,U3,Pb )
%CONSERVEDTOFLUX Summary of this function goes here
%   Detailed explanation goes here

gam = 1.4;


[rho,u,P] = toprimitivevars(U1,U2,U3);
P = Pb;


% E = U3/U1; 



% H = E+P/rho;

% if((U1==0) && (U2 ==0))
%    H = 0; 
%    rho = 0;
%    u = 0;
% end

% % 
% % 
% % J = [ 0 1 0;
% %     0.5*(gam-3)*u^2 (3-gam)*u gam-1;
% %     u*(0.5*(gam-1)*u^2-H) H-(gam-1)*u^2 gam*u];
% % 
% % 
% % 
% % 
% % 
% % F = J*[U1; U2; U3];
% % F1 = F(1);
% % F2 = F(2);
% % F3 = F(3);
E = P/(rho*(gam-1))+0.5*u^2;

F1 = rho*u;
F2 = rho*u^2+P;
F3 = rho*u*(E+P/rho);


end

