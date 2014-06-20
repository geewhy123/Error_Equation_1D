function [ rho,u,P ] = toprimitivevars( U1,U2,U3)
%TOPRIMITIVEVARS Summary of this function goes here
%   Detailed explanation goes here

gam = 1.4;
rho = U1;
u = U2/U1;
P = (gam-1)*(U3-0.5*rho*u^2);

if(U1==0 && U2==0)
   u=0;
   rho=0;
end
end

