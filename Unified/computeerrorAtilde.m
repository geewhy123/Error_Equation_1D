function [ barAtilde ] = computeAtilde( Ut1,Ut2,Ut3 )
%COMPUTEATILDE Summary of this function goes here
%   Detailed explanation goes here
gam = 1.4;
u = Ut2/Ut1;
c = sqrt(gam*(gam-1)*(Ut3/Ut1-0.5*Ut2^2/Ut1^2));
l1 = u;
l2 = u+c;
l3 = u-c;

% [rhot,ut,Pt] = toprimitivevars(Ut1,Ut2,Ut3);
rhot = Ut1;
ut = Ut2;
Pt = Ut3;


ct = sqrt(gam*Pt/rhot);
if(Pt <=0 || rhot <=0)
   ct = 0; 
end




dUdV = eye(3);%[1 0 0; ut rhot 0 ; ut^2/2 rhot*ut 1/(gam-1)];
% dVdU = inv(dUdV);
dVdU = eye(3);%[1 0 0; -ut/rhot 1/rhot 0; (gam-1)*ut^2/2 -(gam-1)*ut (gam-1)];



% Xr = [ut (ut+ct)/(rhot/ct) (ut-ct)/(-rhot/ct); 0 1 1; 0 rhot*ct -rhot*ct];
Xr = [1 rhot/ct -rhot/ct; 0 1 1; 0 rhot*ct -rhot*ct];
Xl = [1 0 -1/ct^2; 0 0.5 1/(2*rhot*ct); 0 0.5 -1/(2*rhot*ct)];



L = abs([l1 0 0; 0 l2 0; 0 0 l3]);

% [l1 l2 l3]
% c
% error('1')

barAtilde = dUdV*Xr*L*Xl*dVdU;



end

