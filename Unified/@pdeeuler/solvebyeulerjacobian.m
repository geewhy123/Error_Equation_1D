function  [errerr2,x,cverr2,exacterr,ee,te  ]  = solvebyeulerjacobian( obj)
%COMPUTEEULERJACOBIAN Summary of this function goes here
%   Detailed explanation goes here
gam = obj.gamma;
ue = obj.exactSolution;
p = obj.pOrder;
q = obj.qOrder;
r = obj.rOrder;
N = obj.nCells;
x = obj.cellCentroids;


V = obj.initialSolution;
U = NaN*ones(3*N+2,1);
obj.computeprimalpseudo();


%%%%%% higher order near bdy
% [Z] = obj.unstructuredrecon(V,p,'solution');

% obj.hOrder = 5;
% m = obj.hOrder ;
% obj.computehigherpseudo();
% 
%       [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
%                [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
%                 [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
%                Zm = [Z1; Z2;Z3];
% 
% error('1')



 %truncation error need exact sol
% J = computeeulerfluxjacobian(obj,ue,'solution');%,x,h,N,p);
% 
% [Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);
% 
% %   [er]=obj.reconplot(Z,'solution')%x,h,N,p,Z);
% %   error('1')
% f = obj.source;
%  [tau]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)
% 
%   tau
% te1 = sum(abs(tau(2:N+1)))/N 
% 
% te = tau;
% % errerr2= NaN;
% % cverr2 = NaN;
% % exacterr = NaN;
% % ee = NaN;
% % return;
% %   error('1')
%  
%  max(abs(tau))
 
 %truncation error need exact sol
  
 del = ones(3*N,1);
 R = ones(3*N+2,1);
 t=0;
 u = ue;
 count = 0;
 c2 = 10;
 dt = 0.0001;

 Rold = R;
 dtold = 1;
 
 while(max(abs(R(2:3*N+1))) > 1e-13 )
     J = obj.computeeulerfluxjacobian(V,'solution');%,x,h,N,p);
    
     count = count +1;
         
     Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 
     dt = dtold*c2*Rratio;

%  spy(J)
%  J
%   error('1')


 K = J(2:3*N+1,2:3*N+1)+eye(3*N)/dt;

[Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
% %%%%%%
% % [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
%   [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
%   [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
%   [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
%                Z = [Z1; Z2;Z3];
% 
% 
%%%%%%

Rold = R;
 [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
 
R(2:3:3*N-1) = phi1(2:N+1);
R(3:3:3*N) = phi2(2:N+1);
R(4:3:3*N+1) = phi3(2:N+1);
u = NaN*ones(N+2,3);
for j = 2:N+1
[u(j,1),u(j,2),u(j,3)] = toconservedvars(V(j,1),V(j,2),V(j,3));
end


    del = K\-R(2:3*N+1);
    
    max(abs(R(2:3*N+1)))

U(2:3:3*N-1) = u(2:N+1,1);
U(3:3:3*N) = u(2:N+1,2);
U(4:3:3*N+1) = u(2:N+1,3);


UU = U(2:3*N+1) + del;
     U = NaN*ones(3*N+2,1);
     U(2:3*N+1) = UU;
     
     
     u(2:N+1,1) = U(2:3:3*N-1) ;
u(2:N+1,2) = U(3:3:3*N) ;
u(2:N+1,3) = U(4:3:3*N+1);

for j = 2:N+1
[V(j,1),V(j,2),V(j,3)] = toprimitivevars(u(j,1),u(j,2),u(j,3));
end
     
     
     t = t+dt;
     
     
     dtold = dt;
     
 end
 
 
%   u
%   error('1')
%   max(abs(u-ue))

[V]
figure
plot(x,V(:,1),'o',x,V(:,2),'v',x,V(:,3),'+')

% % % vv = u-ue;
%   cverr1 = sum(abs(vv(2:N+1)))/N
% % %   cverr2 = sqrt(sum((vv(2:N+1)).^2)/N)

% % %    plot(x,u,'*',x,ue,'o')


 obj.convSoln = u;


Ve = obj.exactSolutionV;
rho = obj.exactSolutionV(:,1);
u = obj.exactSolutionV(:,2);
P = obj.exactSolutionV(:,3);
% error('1')
entropy = log(V(:,3)./(V(:,1).^gam));
   figure
   subplot(2,4,1)
   plot(x,V(:,1),'*',x,rho,'o')
   xlabel('$\rho$','Interpreter','Latex')
      subplot(2,4,2)
   plot(x,V(:,2),'*',x,u,'o')
   xlabel('u')
      subplot(2,4,3)
   plot(x,V(:,3),'*',x,P,'o')
 xlabel('P')
      subplot(2,4,4)
   plot(x,entropy,'*')
xlabel('entropy')

subplot(2,4,5)
plot(x,rho-V(:,1),'x')
subplot(2,4,6)
plot(x,u-V(:,2),'x')
subplot(2,4,7)
plot(x,P-V(:,3),'x')



figure
% plot(x,unew(:,1),x,unew(:,2),x,unew(:,3))
plot(x,V(:,1),x,V(:,2),x,V(:,3))




errerr2 = NaN;
exacterr = Ve-V;
ee = NaN;
cverr2 = [sqrt(sum((exacterr(2:N+1,1)).^2)/N) sqrt(sum((exacterr(2:N+1,2)).^2)/N) sqrt(sum((exacterr(2:N+1,3)).^2)/N)]


Ue = zeros(N+2,3);
for i = 2:N+1
    [Ue(i,1),Ue(i,2),Ue(i,3)] = toconservedvars(Ve(i,1),Ve(i,2),Ve(i,3));
end
exacterr = Ue-obj.convSoln;
cverru2 = [sqrt(sum((exacterr(2:N+1,1)).^2)/N) sqrt(sum((exacterr(2:N+1,2)).^2)/N) sqrt(sum((exacterr(2:N+1,3)).^2)/N)]



te = NaN;





 u=obj.convSoln;
obj.convVreconp = Z;

% Z
% obj.convVreconp
% error('1')




% error('1')


obj.computeprimalleftright();



% obj.convVleft
% obj.convVright
% error('1')


%%%%residual

if(q> 0 && r>0)

  obj.computerespseudo();
  
 

  
  
  Zr = obj.unstructuredrecon(V,r,'residual');

figure
obj.reconplot(Zr(1:r,:),'residual')
figure
obj.reconplot(Zr(r+1:2*r,:),'residual')
figure
obj.reconplot(Zr(2*r+1:3*r,:),'residual')
%

% [Zr] = obj.unstructuredrecon(u,r,'residual');
%   figure
%   obj.reconplot(Zr,'residual');
% error('1')
  
  
  
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
[R1, R2, R3] = obj.computeeulerfluxintegral(Zr,'residual');
 
[R1 R2 R3]

norm1R  = [sum(abs(R1(2:N+1)))/N sum(abs(R2(2:N+1)))/N sum(abs(R3(2:N+1)))/N]

    


 


obj.computeerrorpseudo();



 Zu = obj.unstructuredrecon(V,obj.qOrder,'error');
 obj.convSolnRecon = Zu;
 
 
 
 
%%%error equation

if(obj.bcLeftType == 'D')
   obj.T0 = 0; 
   obj.P0 = 0;
end
if(obj.bcRightType == 'D')
    obj.Pb = 0;
end

exacterrv = obj.exactSolutionV-V
exacterru = obj.exactSolutionU-u
% obj.exactSolutionU
% error('1')

% figure
% plot(x,exacterr,'x')
%  error('1')

% % % need exact stuff
% % % %  obj.computeerrorpseudo();
% % % [Z] = obj.unstructuredrecon(exacterr,q,'error');%ue,x,h,N,NaN,NaN,p);
% % % 
% % % % Z
% % % % figure
% % % % obj.reconplot(Z,'error')
% % % % hold on
% % % % plot(x,exacterr)
% % % % error('1')

f = -[R1 R2 R3];
   obj.errorSource = f;%tau;
 
 
% % %  Je = obj.computefluxjacobian(exacterr,'error');
% % % % obj.errorRM
% % % % error('1')
% % % 
% % % 
% % % 
% % % % error('1')
% % % %  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,obj)
% % % %  error('1')
% % %  [tauE]= obj.computefluxintegral(Z,'error')
% % % %   error('1')
 


 del = ones(3*N,1);
 R = ones(3*N+2,1);
 t=0;
 
 e = exacterrv;
 ee = ones(3*N+2,1);
 ee(1)=NaN;
 ee(3*N+2) = NaN;
 count = 0;
 E = NaN*ones(3*N+2,1);
  dt = 0.001;
dtold = dt;
  c2 = 1;
%  if(obj.qOrder > 4)
% % error('1')
%     kk = 0.00005;
%  end

%  e = 1e-3*ones(size(e));

% for i = 2:N+1
% e(i) = e(i) +1e-4*sin(2*pi*x(i));
% end

%   obj.errorSource = obj.errorSource*0;
%   e = e*0;

% error('1')
%use 0 source, still get NaNs...
% e
% error('1')
 while(max(abs(R(2:3*N+1))) > 1e-13 && count < 10 )
     
    Je = obj.computeeulerfluxjacobian(e,'error');%,x,h,N,p);
%      Jue = obj.computeeulerfluxjacobian(u+e,'error');%,x,h,N,p);
%      Ju = obj.computeeulerfluxjacobian(u,'error');%,x,h,N,p);
% Jue 
% Ju
 Je
%   spy(Je)
%      error('1')
%      Je = Jue-Ju
% e
%      Je
       error('1')
     count = count +1;
%      if(count < 50)
%         dt = kk*(40/N)^2; 
        
      Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 

           dt = 0.001;%dtold*c2*Rratio;
% Rratio


 K = Je(2:3*N+1,2:3*N+1)+eye(3*N)/dt;


[Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
Rold = R;
 [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
    
 R(2:3:3*N-1) = phi1(2:N+1);
R(3:3:3*N) = phi2(2:N+1);
R(4:3:3*N+1) = phi3(2:N+1);
eu = NaN*ones(size(3*N+2,1));
for j = 2:N+1
    eu(j,1) = e(j,1);
    eu(j,2) = e(j,2);
    eu(j,3) = e(j,3);
% [eu(j,1),eu(j,2),eu(j,3)] = toconservedvars(e(j,1),e(j,2),e(j,3));
end

 K
if(isnan(K) )%|| (norm(K) > 1e3))
e
count
Rratio
error('1')
end
 Je
 K
% dt
% Rratio
% dtold
% R

%   error('1')
 del = K\-R(2:3*N+1)

%   error('1')  
%     if(mod(count,100)==0)
    max(abs(R(2:3*N+1)))
%     end
    
E(2:3:3*N-1) = eu(2:N+1,1);
E(3:3:3*N) = eu(2:N+1,2);
E(4:3:3*N+1) = eu(2:N+1,3);


     EE = E(2:3*N+1) + del;%*dt;
     E = NaN*ones(N+2,1);
     E(2:3*N+1) = EE;
     
     eu(2:N+1,1) = E(2:3:3*N-1) ;
     eu(2:N+1,2) = E(3:3:3*N) ;
     eu(2:N+1,3) = E(4:3:3*N+1);

for j = 2:N+1
    e(j,1) = eu(j,1);
    e(j,2) = eu(j,2);
    e(j,3) = eu(j,3);
% [e(j,1),e(j,2),e(j,3)] = toprimitivevars(eu(j,1),eu(j,2),eu(j,3));
end
     
     e
exacterr
del
% error('1')
     
     t = t+dt;
     dtold = dt;
 end

ee = e;
ee
size(exacterr)
size(ee)
w = exacterr-ee

% % % %  Je
% % % %  error('1')
% % % w = Je(2:N+1,2:N+1)\tauE(2:N+1)
% % % 
% % % % x1 = ones(N,1);
% % % %  null(Je(2:N+1,2:N+1),'r')
% % % % Je
% % % %   error('1')
% % % if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(Je(2:N+1,2:N+1)))) < 1e-5)
% % % %   error('5')
% % % %   [Q,R] = qr(Je(2:N+1,2:N+1)') ;
% % % % w = Q*(R'\tauE(2:N+1)) ;
% % % % w = Je(2:N+1,2:N+1)'*Je(2:N+1,2:N+1)\Je(2:N+1,2:N+1)'*tauE(2:N+1)
% % % 
% % % % w = w-dot(w,x1')*x1*dot(w,w);
% % %     w = pinv(Je(2:N+1,2:N+1))*tauE(2:N+1)
% % %   end


max(abs(w))
% error('1')
figure
plot(x,w)

% error('2')
 

% w(2:N+1) = w;
% w(1) = NaN;
% w(N+2) = NaN;
% cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);

ee = exacterr - w;
errerr2 = sqrt(sum((exacterr(2:N+1)-ee(2:N+1)).^2)/N) 
w
figure
subplot(3,1,1)
plot(x,ee(:,1),'*',x,exacterr(:,1),'o')
subplot(3,1,2)
plot(x,ee(:,2),'*',x,exacterr(:,2),'o')
subplot(3,1,3)
plot(x,ee(:,3),'*',x,exacterr(:,3),'o')

else
   errerr2 = NaN;
   exacterr = NaN;
   ee = NaN;
    
end








end

