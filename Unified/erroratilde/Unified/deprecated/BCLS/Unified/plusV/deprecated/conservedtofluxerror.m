function [ F1,F2,F3 ] = conservedtofluxerror( U1,U2,U3 )
%CONSERVEDTOFLUX Summary of this function goes here
%   Detailed explanation goes here
error('1')
gam = 1.4;


[rho,u,P] = toprimitivevars(U1,U2,U3);
E = U3/U1; 



H = E+P/rho;

if((U1==0) || (U2 ==0))
   H = 0; 
   rho = 0;
   u = 0;
end



J = [ 0 1 0;
    0.5*(gam-3)*u^2 (3-gam)*u gam-1;
    u*(0.5*(gam-1)*u^2-H) H-(gam-1)*u^2 gam*u];

% J = [0 1 0;
% 0.5*(gam-3)*U2^2/U1^2 (3-gam)*U2/U1 (gam-1);
% -gam*U2*U3/U1^2+(gam-1)*U2^3/U1^3 gam*U3/U1-1.5*(gam-1)*U2^2/U1^2 gam*U2/U1];
    
% if(abs(U1+0.0949)<1e-3)
% J
% u
% U2/U1
% error('1')
% end
F = J*[U1; U2; U3];
F1 = F(1);
F2 = F(2);
F3 = F(3);



end

