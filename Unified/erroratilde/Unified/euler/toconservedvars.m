function [ u1,u2,u3 ] = toconservedvars( rho,u,P )
%TOCONSERVEDVARS Summary of this function goes here
%   Detailed explanation goes here



u1 = rho;

u2 = rho.*u;
gam = 1.4;
E = (1/(gam-1))*(P./rho)+0.5*u.^2;


if((u==0) && (rho ==0))
   E = 0; 
end


u3 = rho.*E;



% u1 = rho;
% u2 = u;
% u3 = P;



end

