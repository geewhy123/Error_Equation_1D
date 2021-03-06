function [ phi1,phi2,phi3 ] = computeeulerdifffluxintegral( obj,Z,eqn )
%COMPUTEEULERFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here



N=obj.nCells;
h = obj.cellWidths;
x = obj.cellCentroids;
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
rhob= zeros(N+2,1);
ub= zeros(N+2,1);
Pb= zeros(N+2,1);


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




% error('1')
% U = obj.convSoln;
% 
% for i = 2:N+1
% % [V1l(i),V2l(i),V3l(i)]=toprimitivevars(U1l(i),U2l(i),U3l(i));
% % [V1r(i),V2r(i),V3r(i)]=toprimitivevars(U1r(i),U2r(i),U3r(i));
% [V(i,1),V(i,2),V(i,3)]=toprimitivevars(U(i,1),U(i,2),U(i,3));
% V(i,1) = V(i,1) + 0.5*(U1l(i)+U1r(i));
% V(i,2) = V(i,2) + 0.5*(U2l(i)+U2r(i));
% V(i,3) = V(i,3) + 0.5*(U3l(i)+U3r(i));
% 
% end
% fprintf('hack avg');
% % V
% % error('1')
%  Z = obj.unstructuredrecon(V,order,eqn);
% 
% 
% 
% U1l = U1l + obj.convUleft(:,1);
% U2l = U2l + obj.convUleft(:,2);
% U3l = U3l + obj.convUleft(:,3);
% U1r = U1r + obj.convUright(:,1);
% U2r = U2r + obj.convUright(:,2);
% U3r = U3r + obj.convUright(:,3);
% 



for i = 2:N+1
    
    for k = 1:order
    rhor(i) = rhor(i)+ Z(k,i)*(h(i)/2)^(k-1);
    rhol(i) = rhol(i)+ Z(k,i)*(-h(i)/2)^(k-1);
    ur(i)   = ur(i)+ Z(k+order,i)*(h(i)/2)^(k-1);
    ul(i)   = ul(i)+ Z(k+order,i)*(-h(i)/2)^(k-1);
    Pr(i)   = Pr(i)+ Z(k+2*order,i)*(h(i)/2)^(k-1);
    Pl(i)   = Pl(i)+ Z(k+2*order,i)*(-h(i)/2)^(k-1);


    rhob(i) = rhob(i) + (Z(k,i))*obj.moments(i,k);
    ub(i) = ub(i) + (Z(k+order,i))*obj.moments(i,k);
    Pb(i) = Pb(i) + (Z(k+2*order,i))*obj.moments(i,k);
    end

    cl(i) = sqrt(gam*Pl(i)/rhol(i));
    cr(i) = sqrt(gam*Pr(i)/rhor(i));
end

% error('1')



rhol = rhol + obj.convVleft(:,1);
ul = ul + obj.convVleft(:,2);
Pl = Pl + obj.convVleft(:,3);
rhor = rhor + obj.convVright(:,1);
ur = ur + obj.convVright(:,2);
Pr = Pr + obj.convVright(:,3);



 if(strcmp(eqn,'error')==1)
 Vpe = [rhol rhor ul ur Pl Pr];
 V = [obj.convVleft(:,1) obj.convVright(:,1) obj.convVleft(:,2) obj.convVright(:,2) obj.convVleft(:,3) obj.convVright(:,3)];
for i = 2:N+1
 [U(i,1),U(i,2),U(i,3)] = toconservedvars(V(i,1),V(i,3),V(i,5));
end
U;
% error('1')
 end



% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = obj.primalT0; 
% % % %    obj.P0 = obj.primalP0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = obj.primalPb;
% % % % end
U = obj.convSoln;
V = zeros(N+2,3);
for i = 2:N+1
[V(i,1),V(i,2),V(i,3)]=toprimitivevars(U(i,1),U(i,2),U(i,3));
end
V(:,1) = V(:,1) + rhob;
V(:,2) = V(:,2) + ub;
V(:,3) = V(:,3) + Pb;
% V
% error('1')

Z = obj.unstructuredrecon(V,order,eqn);
% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = 0; 
% % % %    obj.P0 = 0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = 0;
% % % % end



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
% U = [U1l U1r U2l U2r U3l U3r]


%%%
% error('1')
% U = obj.convSoln;
% 
% for i = 2:N+1
% % [V1l(i),V2l(i),V3l(i)]=toprimitivevars(U1l(i),U2l(i),U3l(i));
% % [V1r(i),V2r(i),V3r(i)]=toprimitivevars(U1r(i),U2r(i),U3r(i));
% [V(i,1),V(i,2),V(i,3)]=toprimitivevars(U(i,1),U(i,2),U(i,3));
% V(i,1) = V(i,1) + 0.5*(U1l(i)+U1r(i));
% V(i,2) = V(i,2) + 0.5*(U2l(i)+U2r(i));
% V(i,3) = V(i,3) + 0.5*(U3l(i)+U3r(i));
% 
% end
% fprintf('hack avg');
% % V
% % error('1')
% Z = obj.unstructuredrecon(V,order,eqn);
% 
% 
% 
% U1l = U1l + obj.convUleft(:,1);
% U2l = U2l + obj.convUleft(:,2);
% U3l = U3l + obj.convUleft(:,3);
% U1r = U1r + obj.convUright(:,1);
% U2r = U2r + obj.convUright(:,2);
% U3r = U3r + obj.convUright(:,3);

%%%
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
%  [F1l F1r F2l F2r F3l F3r]
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

% l1 = (ur(i)+ul(i+1))/2;
% l2 = (ur(i)+cr(i)+ul(i+1)+cl(i+1))/2;
% l3 = (ur(i)-cr(i)+ul(i+1)-cl(i+1))/2;

% % u = Ut2(i)/Ut1(i);
% % c = sqrt(gam*(gam-1)*(Ut3(i)/Ut1(i)-0.5*Ut2(i)^2/Ut1(i)^2));

% % l1 = u;
% % l2 = u+c;
% % l3 = u-c;
% % 
% % 
% % 
% % if(i==N+1)
% %     u = Ut2(i)/Ut1(i);
% % c = sqrt(gam*(gam-1)*(Ut3(i)/Ut1(i)-0.5*Ut2(i)^2/Ut1(i)^2));
% % 
% % l1 = u;
% % l2 = u+c;
% % l3 = u-c;
% % % l1 = (ur(i));
% % % l2 = (ur(i)+cr(i));
% % % l3 = (ur(i)-cr(i));
% % end


% L = [l1 l2 l3];
% error('2')
bAtilde(:,:,i) = computeAtilde(Ut1(i),Ut2(i),Ut3(i));
end


%  U1L = U1l(3:6)  
%  U1R = U1r(2:5)
%   U2L = U2l(3:6) 
%  U2R = U2r(2:5)
%   U3L = U3l(3:6)
%  U3R = U3r(2:5)
%  
%  [Ut1 Ut2 Ut3]
 
%  bAtilde
%  error('1')

% % % [U1l(3) U1r(2)]
% % % bAtilde(:,:,2)
% error('2')

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
% Ap = @(x) (50/9)*(Ae-At)*(x-2/5);

A = zeros(N+2,1);
phia1 = zeros(N+2,1);
phia2 = zeros(N+2,1);
phia3 = zeros(N+2,1);

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
%     A(i) = (25/9)*(Ae-At)*(xr-2/5)^2+At; 
end
A(1) = obj.getArea(0);
% A(1) = (25/9)*(Ae-At)*(-2/5)^2+At; 



% A
% error('1')
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
   
 phia1(i) =(A(i)*FrAve(i,1)-A(i-1)*FlAve(i,1))/h(i);
phia2(i) = (A(i)*FrAve(i,2)-A(i-1)*FlAve(i,2))/h(i)- PAp(i);
 phia3(i) =(A(i)*FrAve(i,3)-A(i-1)*FlAve(i,3))/h(i);
end











U1l = obj.convUleft(:,1);
U2l = obj.convUleft(:,2);
U3l = obj.convUleft(:,3);
U1r = obj.convUright(:,1);
U2r = obj.convUright(:,2);
U3r = obj.convUright(:,3);
U = obj.convSoln;
for i = 2:N+1
[V1l(i),V2l(i),V3l(i)]=toprimitivevars(U1l(i),U2l(i),U3l(i));
[V1r(i),V2r(i),V3r(i)]=toprimitivevars(U1r(i),U2r(i),U3r(i));
[V(i,1),V(i,2),V(i,3)]=toprimitivevars(U(i,1),U(i,2),U(i,3));
end




% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = obj.primalT0; 
% % % %    obj.P0 = obj.primalP0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = obj.primalPb;
% % % % end

Z = obj.unstructuredrecon(V,order,eqn);

% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = 0; 
% % % %    obj.P0 = 0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = 0;
% % % % end




% VV = [V1l' V1r' V2l' V2r' V3l' V3r']
% UU = [U1l U1r U2l U2r U3l U3r]
% error('1')


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


 if(strcmp(eqn,'error')==1 )%&& obj.T0 == 0)
  [F1l F1r F2l F2r F3l F3r];
%  U
%  [rhol rhor ul ur Pl Pr]
%    error('1')
 end
 
 
 %faces
% [ur ul cr cl]
% error('2')
for i = 2:N
[Ut1(i),Ut2(i),Ut3(i)]=computeroeavg(U1l(i+1),U2l(i+1),U3l(i+1),U1r(i),U2r(i),U3r(i));    

% [Ut1 Ut2 Ut3]


% end check

% l1 = (ur(i)+ul(i+1))/2;
% l2 = (ur(i)+cr(i)+ul(i+1)+cl(i+1))/2;
% l3 = (ur(i)-cr(i)+ul(i+1)-cl(i+1))/2;

% % u = Ut2(i)/Ut1(i);
% % c = sqrt(gam*(gam-1)*(Ut3(i)/Ut1(i)-0.5*Ut2(i)^2/Ut1(i)^2));

% % l1 = u;
% % l2 = u+c;
% % l3 = u-c;
% % 
% % 
% % 
% % if(i==N+1)
% %     u = Ut2(i)/Ut1(i);
% % c = sqrt(gam*(gam-1)*(Ut3(i)/Ut1(i)-0.5*Ut2(i)^2/Ut1(i)^2));
% % 
% % l1 = u;
% % l2 = u+c;
% % l3 = u-c;
% % % l1 = (ur(i));
% % % l2 = (ur(i)+cr(i));
% % % l3 = (ur(i)-cr(i));
% % end


% L = [l1 l2 l3];
% error('2')
bAtilde(:,:,i) = computeAtilde(Ut1(i),Ut2(i),Ut3(i));
end


%  U1L = U1l(3:6)  
%  U1R = U1r(2:5)
%   U2L = U2l(3:6) 
%  U2R = U2r(2:5)
%   U3L = U3l(3:6)
%  U3R = U3r(2:5)
%  
%  [Ut1 Ut2 Ut3]
 
%  bAtilde
%  error('1')

% % % [U1l(3) U1r(2)]
% % % bAtilde(:,:,2)
% error('2')

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










% error('2')
PAp = NaN*ones(N+2,1);
% Ap = @(x) (50/9)*(Ae-At)*(x-2/5);

A = zeros(N+2,1);
phib1 = zeros(N+2,1);
phib2 = zeros(N+2,1);
phib3 = zeros(N+2,1);
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
%     A(i) = (25/9)*(Ae-At)*(xr-2/5)^2+At; 
end
A(1) = obj.getArea(0);
% A(1) = (25/9)*(Ae-At)*(-2/5)^2+At; 



% A
% error('1')
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
   
 phib1(i) =(A(i)*FrAve(i,1)-A(i-1)*FlAve(i,1))/h(i);
phib2(i) = (A(i)*FrAve(i,2)-A(i-1)*FlAve(i,2))/h(i)- PAp(i);
 phib3(i) =(A(i)*FrAve(i,3)-A(i-1)*FlAve(i,3))/h(i);
end


phi1 = phia1-phib1;
phi2 = phia2-phib2;
phi3 = phia3-phib3;


    phi1 = phi1 - obj.errorSource(:,1);
    phi2 = phi2 - obj.errorSource(:,2);
    phi3 = phi3 - obj.errorSource(:,3);



PP = [phia1 phia2 phia3 phib1 phib2 phib3];

%  error('1')




% [phi1(N+1) phi2(N+1) phi3(N+1)]
% % %  [FlAve(:,:) FrAve(:,:)]

% error('1')


% i = N+1;
% [A(N+1) A(N) FrAve(N+1,1) FrAve(N+1,2) FrAve(N+1,3) FlAve(N+1,1) FlAve(N+1,2) FlAve(N+1,3) ] 
% 0.5*[(F1l(i)+F1r(i-1)); (F2l(i)+F2r(i-1)); (F3l(i)+F3r(i-1))]
% -0.5*bAtilde(:,:,i-1)*([ U1l(i); U2l(i); U3l(i)]- [ U1r(i-1); U2r(i-1); U3r(i-1)]) 


% [Ut1(N+1) Ut2(N+1) Ut3(N+1)]
% error('3')

% El = (1/(gam-1))*Pl./rhol+0.5*ul.^2;
% f1l = rhol.*ul;
% f2l = rhol.*ul.^2+Pl;
% f3l = rhol.*ul.*(El+Pl./rhol);
% [f1l f2l f3l]

% error('1')

% % if(strcmp(eqn,'error')==1)
%     phi1 = phi1 - obj.errorSource(:,1);
%     phi2 = phi2 - obj.errorSource(:,2);
%     phi3 = phi3 - obj.errorSource(:,3);
% end

end

