function [ phi1,phi2,phi3 ] = computeeulerfluxintegral( obj,Z,eqn )
%COMPUTEEULERFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
N=obj.nCells;
h = obj.cellWidths;
rhol  = zeros(N+2,1);
rhor  = zeros(N+2,1);
ul = zeros(N+2,1);
ur = zeros(N+2,1);
Pl = zeros(N+2,1);
Pr = zeros(N+2,1);
cr = zeros(N+2,1);
cl = zeros(N+2,1);
Ut1 = zeros(N+2,1);
Ut2 = zeros(N+2,1);
Ut3 =zeros(N+2,1);
gam = 1.4;
bAtilde = zeros(3,3,N+2);
FrAve = zeros(N+2,3);
FlAve = zeros(N+2,3);
Z
for i = 2:N+1
    
    for k = 1:obj.pOrder
    rhor(i) = rhor(i)+ Z(k,i)*(h(i)/2)^(k-1);
    rhol(i) = rhol(i)+ Z(k,i)*(-h(i)/2)^(k-1);
    ur(i)   = ur(i)+ Z(k+obj.pOrder,i)*(h(i)/2)^(k-1);
    ul(i)   = ul(i)+ Z(k+obj.pOrder,i)*(-h(i)/2)^(k-1);
    Pr(i)   = Pr(i)+ Z(k+2*obj.pOrder,i)*(h(i)/2)^(k-1);
    Pl(i)   = Pl(i)+ Z(k+2*obj.pOrder,i)*(-h(i)/2)^(k-1);
    end
    
    cl(i) = sqrt(gam*Pl(i)/rhol(i));
    cr(i) = sqrt(gam*Pr(i)/rhor(i));
end
[rhol rhor ul ur Pl Pr]
% error('1')
[U1l,U2l,U3l]=toconservedvars(rhol,ul,Pl);
[U1r,U2r,U3r]=toconservedvars(rhor,ur,Pr);

[U1l U1r U2l U2r U3l U3r]


F1l = zeros(N+2,1);
F2l = zeros(N+2,1);
F3l = zeros(N+2,1);
F1r = zeros(N+2,1);
F2r = zeros(N+2,1);
F3r = zeros(N+2,1);

for i = 2:N+1
[F1l(i),F2l(i),F3l(i)]=conservedtoflux(U1l(i),U2l(i),U3l(i));
[F1r(i),F2r(i),F3r(i)]=conservedtoflux(U1r(i),U2r(i),U3r(i));
end

[F1l F1r F2l F2r F3l F3r]
% error('1')
%faces
for i = 2:N+1
[Ut1(i),Ut2(i),Ut3(i)]=computeroeavg(U1l(i+1),U2l(i+1),U3l(i+1),U1r(i),U2r(i),U3r(i));    

% [Ut1 Ut2 Ut3]


% end check

l1 = (ur(i)+ul(i+1))/2;
l2 = (ur(i)+cr(i)+ul(i+1)+cl(i+1))/2;
l3 = (ur(i)-cr(i)+ul(i+1)-cl(i+1))/2;


[l1 l2 l3]
% error('2')
bAtilde(:,:,i) = computeAtilde(Ut1(i),Ut2(i),Ut3(i),l1,l2,l3);
end


for i = 2:N+1
FrAve(i,1:3) =( 0.5*[(F1r(i)+F1l(i+1)); (F2r(i)+F2l(i+1)); (F3r(i)+F3l(i+1))]  -0.5*bAtilde(:,:,i)*  ([ U1l(i+1); U2l(i+1); U3l(i+1)]- [ U1r(i); U2r(i); U3r(i)]) )';
FlAve(i,1:3) =( 0.5*[(F1l(i)+F1r(i-1)); (F2l(i)+F2r(i-1)); (F3l(i)+F3r(i-1))]  -0.5*bAtilde(:,:,i-1)*([ U1l(i); U2l(i); U3l(i)]- [ U1r(i-1); U2r(i-1); U3r(i-1)]) )';
end

FlAve(2,1:3) = [F1l(2); F2l(2); F3l(2)]';
% FrAve(1:3,2) = [0;0;0];
% FlAve(1:3,N+1)= [ 0;0;0];
FrAve(N+1,1:3) = [F1r(N+1);F2r(N+1);F3r(N+1)]';

[FlAve FrAve]
error('2')
A = ones(N+2,1);
P = A;
Ap = A;

for i = 2:N+1
   
 phi1(i) =(A(i)*FrAve(i,1)-A(i-1)*FlAve(i,1))/h(i);
phi2(i) = (A(i)*FrAve(i,2)-A(i-1)*FlAve(i,2))/h(i)- P(i)*Ap(i);
 phi3(i) =(A(i)*FrAve(i,3)-A(i-1)*FlAve(i,3))/h(i);
end

% El = (1/(gam-1))*Pl./rhol+0.5*ul.^2;
% f1l = rhol.*ul;
% f2l = rhol.*ul.^2+Pl;
% f3l = rhol.*ul.*(El+Pl./rhol);
% [f1l f2l f3l]

% error('1')
end

