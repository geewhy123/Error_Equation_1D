function [ F1,F2,F3 ] = outboundaryfluxNewton( U1,U2,U3,e1,e2,e3,Pb )
%CONSERVEDTOFLUX Summary of this function goes here
%   Detailed explanation goes here

gam = 1.4;





[rho,u,P] = toprimitivevars(U1,U2,U3);
% P = Pb;


dVdU = [1 0 0;
    -u/rho 1/rho 0;
    ((gam-1)/2)*u^2 -(gam-1)*u gam-1];
dUdV = inv(dVdU);

enew = dVdU*[e1 e2 e3]';

ev1 = enew(1);
ev2 = enew(2);
ev3 = enew(3);

ev3 = 0;

eunew = dUdV*[ev1 ev2 ev3]';
e1 = eunew(1);
e2 = eunew(2);
e3 = eunew(3);
Q = zeros(3,3);
Q(1,1) = 0;
        Q(1,2) = 1;
        Q(1,3) = 0;
        Q(2,1) = ((gam-1)/2-1)*U2^2/U1^2;
        Q(2,2) = (3-gam)*U2/U1;
        Q(2,3) = (gam-1);
        Q(3,1) = (gam-1)*U2^3/U1^3-gam*U2*U3/U1^2;
        Q(3,2) = gam*U3/U1 - 1.5*(gam-1)*U2^2/U1^2;
        Q(3,3) = gam*U2/U1;
        
      
        
        F = Q*[e1 e2 e3]';
        F1 = F(1);
        F2 = F(2);
        F3 = F(3);
    


return;
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
% [F1 F2 F3]

end

