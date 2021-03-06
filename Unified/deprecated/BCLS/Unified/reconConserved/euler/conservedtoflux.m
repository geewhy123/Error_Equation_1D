function [ F1,F2,F3 ] = conservedtoflux( U1,U2,U3 )
%CONSERVEDTOFLUX Summary of this function goes here
%   Detailed explanation goes here

gam = 1.4;


[rho,u,P] = toprimitivevars(U1,U2,U3);
E = U3/U1; 



H = E+P/rho;

% if((U1==0) && (U2 ==0))
%    H = 0; 
%    rho = 0;
%    u = 0;
% end



J = [ 0 1 0;
    0.5*(gam-3)*u^2 (3-gam)*u gam-1;
    u*(0.5*(gam-1)*u^2-H) H-(gam-1)*u^2 gam*u];





F = J*[U1; U2; U3];
F1 = F(1);
F2 = F(2);
F3 = F(3);


F1 = U2;
F2 = ((3-gam)/2)*(U2^2/U1)+(gam-1)*U3;
F3 = gam*U2*U3/U1-0.5*(gam-1)*U2^2/U1^2;
end

