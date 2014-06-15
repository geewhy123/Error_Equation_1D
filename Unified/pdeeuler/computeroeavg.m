function [ Ut1,Ut2,Ut3 ] = computeroeavg( Ul1,Ul2,Ul3,Ur1,Ur2,Ur3 )
%COMPUTEROEAVG Summary of this function goes here
%   Detailed explanation goes here

% [Ul1 Ul2 Ul3 Ur1 Ur2 Ur3]
% error('3')
gam = 1.4;
[rhor, ur,Pr]= toprimitivevars(Ur1,Ur2,Ur3);
Er = (1/(gam-1))*(Pr./rhor)+0.5*ur.^2;
hr = Er+Pr/rhor;

[rhol, ul,Pl]= toprimitivevars(Ul1,Ul2,Ul3);
El = (1/(gam-1))*(Pl./rhol)+0.5*ul.^2;
hl = El+Pl/rhol;

rhot = sqrt(rhor*rhol);

ut = (sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol));

ht = (sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol));


Pt =(ht-0.5*ut^2)*((gam-1)/(gam))*rhot;


Ut1 = rhot;
Ut2 = rhot*ut;
Ut3 = rhot*ht-Pt;



end

