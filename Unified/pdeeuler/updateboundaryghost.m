function [ output_args ] = updateboundaryghost( U,h,k,N)
%UPDATEBOUNDARYGHOST Summary of this function goes here
%   Detailed explanation goes here

gam = 1.4;

rho1 = U(2,1);
u1 = U(2,2);
P1 = U(2,3);

ug = U(1,2);
Pg = U(1,3);

c1 = sqrt(gam*P1/rho1);

l3 = (u1-c1)*k/h;
P0 = 1;
T0 = 1;

%cv = 1;
cv = 1/(gam-1);
cs = 2*gam*((gam-1)/(gam+1))*cv*T0;
dPdu = P0*(gam/(gam-1))*(1-((gam-1)/(gam+1))*u1^2/cs^2)^(1/(gam-1))*(-2*(u1/cs)*(gam-1)/(gam+1));



% Pg = P0*(1-((gam-1)/(gam+1))*ug^2/cs^2)^(gam/(gam-1));



deltaug = (-l3/(1-l3))*(P1-Pg-rho1*c1*(u1-ug))/(dPdu-rho1*c1);

ugnew = ug+deltaug;

Pgnew = P0*(1-((gam-1)/(gam+1))*ugnew^2/cs^2)^(gam/(gam-1));
Tgnew = T0*(1-((gam-1)/(gam+1))*ugnew^2/cs^2);

R = (gam-1)*cv;

rhognew = Pgnew/(R*Tgnew);

U(1,1) = rhognew ;
U(1,2) = ugnew;
U(1,3) = Pgnew;


%outflow
rhoI = U(N+1,1);
uI = U(N+1,2);
PI = U(N+1,3);
rhog = U(N+2,1);
ug = U(N+2,2);
Pg = U(N+2,3);
cI = sqrt(gam*PI/rhoI);

l1 = (uI);
l2 = (uI+cI)*k/h;

rho = rhoI;
c = cI;

deltaPg = 0;
deltarhog = (1/c^2)*deltaPg-(l1/(1+l1))*(rhog-rhoI-(1/c^2)*(Pg-PI));
deltaug = (1/(rho*c))*(-deltaPg-(l2/(1+l2))*(Pg-PI+rho*c*(ug-uI)) ) ;

U(N+2,1) = deltarhog;
U(N+2,2) = deltaug;
U(N+2,3) = deltaPg;



end

