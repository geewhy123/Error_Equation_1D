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
gam = 1.4;
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
    
    
end
% [rhol rhor ul ur Pl Pr]
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


% El = (1/(gam-1))*Pl./rhol+0.5*ul.^2;
% f1l = rhol.*ul;
% f2l = rhol.*ul.^2+Pl;
% f3l = rhol.*ul.*(El+Pl./rhol);
% [f1l f2l f3l]

error('1')
end

