function [ phi1,phi2,phi3 ] = computeeulerfluxintegral( obj,Z,eqn )
%COMPUTEEULERFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
if(strcmp(eqn,'error')==1)
%    [ phi1,phi2,phi3 ] = computeeulererrorfluxintegral( obj,Z,eqn )
   [ phi1,phi2,phi3 ] = computeeulerdifffluxintegral( obj,Z,eqn );
    return
end


N=obj.nCells;
h = obj.cellWidths;
x = obj.cellCentroids;
sourceMMS = obj.source;
rhol  = zeros(N+2,1);
rhor  = zeros(N+2,1);
ul = zeros(N+2,1);
ur = zeros(N+2,1);
U1l = zeros(N+2,1);
U2l = zeros(N+2,1);
U3l = zeros(N+2,1);
U1r = zeros(N+2,1);
U2r = zeros(N+2,1);
U3r = zeros(N+2,1);

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

if(strcmp(eqn,'solution')==1)
   order = obj.pOrder; 
elseif(strcmp(eqn,'residual')==1)
   order = obj.rOrder;
elseif(strcmp(eqn,'error')==1)
   order = obj.qOrder;
else
    assert(0)
end


for i = 2:N+1

%%%%higher
if(obj.hOrder > 0)%|| (strcmp(eqn,'residual')==1 && obj.hOrder > obj.rOrder) || (strcmp(eqn,'error')==1 && obj.hOrder > obj.qOrder) )
if(i==2 || i == 3 || i == N || i == N+1)
    order = obj.hOrder;
end
end
% %%%%higher


    for k = 1:order
    rhor(i) = rhor(i)+ Z(k,i)*(h(i)/2)^(k-1);
    rhol(i) = rhol(i)+ Z(k,i)*(-h(i)/2)^(k-1);
    ur(i)   = ur(i)+ Z(k+order,i)*(h(i)/2)^(k-1);
    ul(i)   = ul(i)+ Z(k+order,i)*(-h(i)/2)^(k-1);
    Pr(i)   = Pr(i)+ Z(k+2*order,i)*(h(i)/2)^(k-1);
    Pl(i)   = Pl(i)+ Z(k+2*order,i)*(-h(i)/2)^(k-1);
    end
    
    cl(i) = sqrt(gam*Pl(i)/rhol(i));
    cr(i) = sqrt(gam*Pr(i)/rhor(i));
end

% if(strcmp(eqn,'error')==1)
 V = [rhol rhor ul ur Pl Pr];
% end

% rhor(N+1) = obj.bcRightVal(1);
% ur(N+1) = obj.bcRightVal(2);
% Pr(N+1) = obj.Pb;

%%%

% Z
% error('1')

if(strcmp(eqn,'solution') && (min(rhol)<0 || min(rhor)<0 || min(Pl)<0 || min(Pr)<0))
   figure
obj.reconplot(Z(1:order,:),'solution')

obj.reconplot(Z(order+1:2*order,:),'solution')

obj.reconplot(Z(2*order+1:3*order,:),'solution')
Z
error('2')
   
    
  
 
end
for i = 2:N+1
[U1l(i),U2l(i),U3l(i)]=toconservedvars(rhol(i),ul(i),Pl(i));
[U1r(i),U2r(i),U3r(i)]=toconservedvars(rhor(i),ur(i),Pr(i));
end
U = [U1l U1r U2l U2r U3l U3r];




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


% if(strcmp(eqn,'error')==1 && obj.T0 == 0)
  [F1l F1r F2l F2r F3l F3r];
%  U
%  [rhol rhor ul ur Pl Pr]
%   error('1')
% end
 
 
 %faces
% [ur ul cr cl]
% error('2')
for i = 2:N
[Ut1(i),Ut2(i),Ut3(i)]=computeroeavg(U1l(i+1),U2l(i+1),U3l(i+1),U1r(i),U2r(i),U3r(i));    

% [Ut1 Ut2 Ut3]


% end check


bAtilde(:,:,i) = computeAtilde(Ut1(i),Ut2(i),Ut3(i));
end

save('atilde.mat','bAtilde')


for i = 2:N+1
FrAve(i,1:3) =(0.5*[(F1r(i)+F1l(i+1)); (F2r(i)+F2l(i+1)); (F3r(i)+F3l(i+1))]  -0.5*bAtilde(:,:,i)*  ([ U1l(i+1); U2l(i+1); U3l(i+1)]- [ U1r(i); U2r(i); U3r(i)]) )';
FlAve(i,1:3) =(0.5*[(F1l(i)+F1r(i-1)); (F2l(i)+F2r(i-1)); (F3l(i)+F3r(i-1))]  -0.5*bAtilde(:,:,i-1)*([ U1l(i); U2l(i); U3l(i)]- [ U1r(i-1); U2r(i-1); U3r(i-1)]) )';
end
% % % [FlAve FrAve]
FlAve(2,1:3) = [F1l(2); F2l(2); F3l(2)]';
% FrAve(1:3,2) = [0;0;0];
% FlAve(1:3,N+1)= [ 0;0;0];
FrAve(N+1,1:3) = [F1r(N+1);F2r(N+1);F3r(N+1)]';

F = [FlAve FrAve];
% error('2')
PAp = NaN*ones(N+2,1);


A = zeros(N+2,1);
phi1 = zeros(N+2,1);
phi2 = zeros(N+2,1);
phi3 = zeros(N+2,1);

  c1 = 0.3478548451;
c2 = 0.6521451549;
c3 = 0.6521451549;
c4 = 0.3478548451;
x1= 0.8611363116;
x2 = 0.339981436;
x3 = -0.339981436;
x4= -0.8611363116;


for i = 2:N+1
    xr = x(i)+h(i)/2;   
   A(i)=obj.getArea(xr);
 
end
A(1) = obj.getArea(0);


for i = 2:N+1
    xl = x(i)-h(i)/2;
    xr = x(i)+h(i)/2;
    
    
     xx1 = ((xr-xl)/2)*x1+(xr+xl)/2;
 xx2 = ((xr-xl)/2)*x2+(xr+xl)/2;
 xx3 = ((xr-xl)/2)*x3+(xr+xl)/2;
 xx4 = ((xr-xl)/2)*x4+(xr+xl)/2;
 P1 = 0;
 P2 = 0;
 P3 = 0;
 P4 = 0;
    for k = 1:order
       P1 = P1 + Z(k+2*order,i)*(xx1-x(i))^(k-1); 
       P2 = P2 + Z(k+2*order,i)*(xx2-x(i))^(k-1) ;
       P3 = P3 + Z(k+2*order,i)*(xx3-x(i))^(k-1) ;
       P4 = P4 + Z(k+2*order,i)*(xx4-x(i))^(k-1) ;
    end
 
 
 PAp(i) = (1/h(i))*(c1*P1*obj.getAp(xx1)+c2*P2*obj.getAp(xx2)+c3*P3*obj.getAp(xx3)+c4*P4*obj.getAp(xx4))*(xr-xl)/2;
   
 phi1(i) =(A(i)*FrAve(i,1)-A(i-1)*FlAve(i,1))/h(i) - sourceMMS(i,1);
phi2(i) = (A(i)*FrAve(i,2)-A(i-1)*FlAve(i,2))/h(i)- PAp(i) - sourceMMS(i,2);
 phi3(i) =(A(i)*FrAve(i,3)-A(i-1)*FlAve(i,3))/h(i) - sourceMMS(i,3);
 
end


end

