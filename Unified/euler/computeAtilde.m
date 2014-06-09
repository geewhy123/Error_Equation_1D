function [ barAtilde ] = computeAtilde( Ut1,Ut2,Ut3 ,l1,l2,l3)
%COMPUTEATILDE Summary of this function goes here
%   Detailed explanation goes here
gam = 1.4;
[rhot,ut,Pt] = toprimitivevars(Ut1,Ut2,Ut3);
ct = sqrt(gam*Pt/rhot);

dUdV = [1 0 0; ut rhot 0 ; ut^2/2 rhot*ut 1/(gam-1)];
% dVdU = inv(dUdV);
dVdU = [1 0 0; -ut/rhot 1/rhot 0; (gam-1)*ut^2/2 -(gam-1)*ut (gam-1)];

Xr = [ut (ut+ct)/(rhot/ct) (ut-ct)/(-rhot/ct); 0 1 1; 0 rhot*ct -rhot*ct];
Xl = [1 0 -1/ct^2; 0 0.5 1/(2*rhot*ct); 0 0.5 -1/(2*rhot*ct)];



L = [l1 0 0; 0 l2 0; 0 0 l3];


barAtilde = dUdV*Xr*L*Xl*dVdU;

end

