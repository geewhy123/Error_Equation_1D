function [ Ut1,Ut2,Ut3 ] = computeerrorroeavg(obj, Ul1,Ul2,Ul3,Ur1,Ur2,Ur3,i )
%COMPUTEROEAVG Summary of this function goes here
%   Detailed explanation goes here

N = obj.nCells;
h = obj.cellWidths;
U = zeros(N+2,3);
% primrhor = zeros(N+2,1);
% primur = zeros(N+2,1);
% primPr= zeros(N+2,1);
% primrhol = zeros(N+2,1);
% primul = zeros(N+2,1);
% primPl= zeros(N+2,1);
primrhor = 0;
primrhol = 0;
primur = 0;
primul = 0;
primPl = 0;
primPr = 0;
order = obj.pOrder;
Z = obj.convSolnRecon;

% for i = 2:N+1
    
    for k = 1:order
    primrhor = primrhor+ Z(k,i)*(h(i)/2)^(k-1);
    primrhol = primrhol+ Z(k,i)*(-h(i)/2)^(k-1);
    primur   = primur+ Z(k+order,i)*(h(i)/2)^(k-1);
    primul   = primul+ Z(k+order,i)*(-h(i)/2)^(k-1);
    primPr   = primPr+ Z(k+2*order,i)*(h(i)/2)^(k-1);
    primPl   = primPl+ Z(k+2*order,i)*(-h(i)/2)^(k-1);
    end
    
%     cl(i) = sqrt(gam*Pl(i)/rhol(i));
%     cr(i) = sqrt(gam*Pr(i)/rhor(i));
% end


% [U]=obj.convSoln;

% for i = 2:N+1
% [primrhor(i),primu(i),primP(i)] = toprimitivevars(U(i,1),U(i,2),U(i,3));
% end
%  [primrhol,primrhor,primul,primur,primPl,primPr]
% [primrho,primu,primP] = toprimitivevars(U(:,1),U(:,2),U(:,3));
%  error('5')

%   [Ul1 Ul2 Ul3 Ur1 Ur2 Ur3]
%   error('3')
gam = 1.4;
% [rhor, ur,Pr]= toprimitivevars(Ur1,Ur2,Ur3);
rhor = Ur1;
ur = Ur2;
Pr = Ur3;

Er = (1/(gam-1))*(Pr./rhor)+0.5*ur.^2;
hr = Er+Pr/rhor;
if(rhor == 0 )
    hr = 0;
end


% [rhol, ul,Pl]= toprimitivevars(Ul1,Ul2,Ul3);
rhol = Ul1;
ul = Ul2;
Pl = Ul3;


El = (1/(gam-1))*(Pl./rhol)+0.5*ul.^2;
hl = El+Pl/rhol;

if(rhol == 0)
    hl = 0;
end


% [rhol rhor ul ur Pl Pr]
% error('3')

 rhot = (rhor+rhol)/2;%sqrt(primrhor*primrhol);%sqrt(rhor*rhol);

ut = (ur+ul)/2;%(sqrt(primrhor)*ur+sqrt(primrhol)*ul)/(sqrt(primrhor)+sqrt(primrhol));

ht = (hr+hl)/2;%(sqrt(primrhor)*hr+sqrt(primrhol)*hl)/(sqrt(primrhor)+sqrt(primrhol));


Pt =(ht-0.5*ut^2)*((gam-1)/(gam))*rhot;


Ut1 = rhot;
Ut2 = rhot*ut;
Ut3 = rhot*ht-Pt;


%no conversion...
Ut1 = (Ul1+Ur1)/2;
Ut2 = (Ul2+Ur2)/2;
Ut3 = (Ul3+Ur3)/2;







end

