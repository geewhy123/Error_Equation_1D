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


%%%%% higher order near bdy
% [Z] = obj.unstructuredrecon(V,p,'solution');
% 
obj.hOrder = 0;%obj.pOrder;
% m = obj.hOrder ;
obj.computehigherpseudo();
% 
%       [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
%                [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
%                 [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
%                Zm = [Z1; Z2;Z3];
% 
% % error('1')



 %truncation error need exact sol
teu = zeros(N+2,3);
% tev = zeros(N+2,3);
Ve = obj.exactSolutionV;




% % % % higher
[Z] = obj.unstructuredrecon(Ve,p,'solution');%ue,x,h,N,NaN,NaN,p);
% % % % %  [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
% % % % %                [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
% % % % %                 [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
% % % % %                Z = [Z1; Z2;Z3];
% % % % % % % % % % 



%truncation error

 [tauu1 tauu2 tauu3]=obj.computeeulerfluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

  [tauu1 tauu2 tauu3]
te1 = [ sum(abs(tauu1(2:N+1)))/N  sum(abs(tauu2(2:N+1)))/N sum(abs(tauu3(2:N+1)))/N ]

teu = [tauu1 tauu2 tauu3];
 
%truncation error need exact sol

  
 R = ones(3*N+2,1);
 t=0;
 u = ue;
 count = 0;
 c2 = 10;
 dtold = 0.01;
 Rold = R;
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
% % %%%%%
% % [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
% %   [Z3] = higherunstructuredreconeuler (obj,V(:,3),m,'solution',3);                          
% %   [Z1] = higherunstructuredreconeuler (obj,V(:,1),m,'solution',1); 
% %   [Z2] = higherunstructuredreconeuler (obj,V(:,2),m,'solution',2);
% %                Z = [Z1; Z2;Z3];
% % 
% % 
% % %%%%

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
 
 
figure
plot(x,V(:,1),'o',x,V(:,2),'v',x,V(:,3),'+')


 obj.convSoln = u;


Ve = obj.exactSolutionV;
rho = obj.exactSolutionV(:,1);
u = obj.exactSolutionV(:,2);
P = obj.exactSolutionV(:,3);


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
subplot(2,4,8)
plot(x,V(:,3)./V(:,1).^gam,'+')



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



te = teu;


 u=obj.convSoln;
obj.convSolutionV = V;
 [Z] = obj.unstructuredrecon(V,p,'solution');
obj.convVreconp = Z;


obj.computeprimalleftright();



[phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,'solution')




Pinf = 1;
Tinf = 1;
T0 = obj.T0;
P0 = obj.P0;
rhoinf = Pinf/Tinf;
cinf = sqrt(gam*Pinf/rhoinf);
uinf = cinf*sqrt((2/(gam-1))*(T0/Tinf-1));

rhoi = 0;
ui = 0;
Pi = 0;

h = obj.cellWidths;
for k = 1:p
    
    rhoi = rhoi+ Z(k,2)*(-h(2)/2)^(k-1);
    ui   = ui+ Z(k+p,2)*(-h(2)/2)^(k-1);
    Pi   = Pi+ Z(k+2*p,2)*(-h(2)/2)^(k-1);
end

Vi= [rhoi ui Pi]'
ci = sqrt(gam*Pi/rhoi);

A = [cinf^2 0 -1; 0 -rhoinf*cinf -1; 0 rhoi*ci -1];
b = [rhoinf-Pinf; -rhoinf*cinf*uinf-Pinf; rhoi*ci*ui-Pi];

A\b

A\b-Vi


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



%%%%
[Zq] = obj.unstructuredrecon(Ve,q,'error');%ue,x,h,N,NaN,NaN,p);
obj.reconexactsolutionV = Zq;
%%%%


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

exacterrv = obj.exactSolutionV-V;
exacterru = obj.exactSolutionU-u

obj.exactSolutionV
UU =zeros(N+2,3);
for z = 2:N+1
[UU(z,1),UU(z,2),UU(z,3)] = toconservedvars(obj.exactSolutionV(z,1),obj.exactSolutionV(z,2),obj.exactSolutionV(z,3)); 
end
UU

obj.exactSolutionU
% % error('1')




% error('1')

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


% load('tau4.mat')

f = -[R1 R2 R3];
   obj.errorSource = f;teu;%f;%tau;
   fprintf('using t.e. as source')

   
   figure
   subplot(3,1,1)
   plot(x,f(:,1),'+',x,teu(:,1),'^')
%    figure
   subplot(3,1,2)
   plot(x,f(:,2),'+',x,teu(:,2),'^')
%    figure
subplot(3,1,3)
   plot(x,f(:,3),'+',x,teu(:,3),'^')
   legend('residual source','te source')
   [f teu]
   [ mean(abs(f(2:N+1,1))) mean(abs(f(2:N+1,2))) mean(abs(f(2:N+1,3)))]
%    error('1')
 
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
  
dtold = 0.01;
  c2 = 10;


Upe =zeros(N+2,3);
Vpe = zeros(N+2,3);

%
    [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);
    [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)

%      [Zvpe] = obj.unstructuredrecon(e+obj.convSolutionV,q,'error');%u,x,h,N,NaN,NaN,p);
%     [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Zvpe,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)

    obj.convVleft
    obj.convVright
    
    [phi1 phi2 phi3]
    figure
    plot(x,phi1,x,phi2,x,phi3)
%      error('1')
%
    R(2:3:3*N-1) = phi1(2:N+1);
    R(3:3:3*N) = phi2(2:N+1);
    R(4:3:3*N+1) = phi3(2:N+1);
 while(max(abs(R(2:3*N+1))) > 1e-13  )
     
% % % % if(obj.bcLeftType == 'D')
% % % %    obj.T0 = 0; 
% % % %    obj.P0 = 0;
% % % % end
% % % % if(obj.bcRightType == 'D')
% % % %     obj.Pb = 0;
% % % % end



    Je = obj.computeeulerfluxjacobian(e,'error');%,x,h,N,p);


% % % % % if(obj.bcLeftType == 'D')
% % % % %    obj.T0 = 1; 
% % % % %    obj.P0 = 1;
% % % % % end
% % % % % if(obj.bcRightType == 'D')
% % % % %     obj.Pb = 0.95;
% % % % % end



%      Jue = obj.computeeulerfluxjacobian(u+e,'error');%,x,h,N,p);
%      Ju = obj.computeeulerfluxjacobian(u,'error');%,x,h,N,p);
% Jue 
% Ju
%   Je
%       spy(Je)
%        error('1')
%      Je = Jue-Ju
% e
%      Je
%         error('1')
     count = count +1;
%      if(count < 50)
%         dt = kk*(40/N)^2; 
        
     Rratio =norm(Rold(2:3*N+1),2)/norm(R(2:3*N+1),2); 

     dt= dtold*c2*Rratio;

     K = Je(2:3*N+1,2:3*N+1)+eye(3*N)/dt;




    [Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);


    Rold = R;
    [phi1 phi2 phi3]=obj.computeeulerfluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
    
    R(2:3:3*N-1) = phi1(2:N+1);
    R(3:3:3*N) = phi2(2:N+1);
    R(4:3:3*N+1) = phi3(2:N+1);
    eu = NaN*ones(N+2,3);
%     V
%     obj.convVreconp
%     [obj.convVleft obj.convVright]
% Z
%     [phi1 phi2 phi3]
%      figure
%      plot(x,phi1,x,phi2,x,phi3)
%     figure
%     plot(x,e)
%      error('1')

%%%%%new translation
    Vpe = e+V;


    for j = 2:N+1
%     eu(j,1) = e(j,1);
%     eu(j,2) = e(j,2);
%     eu(j,3) = e(j,3);
% [eu(j,1),eu(j,2),eu(j,3)] = toconservedvars(e(j,1),e(j,2),e(j,3));
        [Upe(j,1),Upe(j,2),Upe(j,3)] = toconservedvars(Vpe(j,1),Vpe(j,2),Vpe(j,3));
    end

    eu = Upe-u;
% error('1')
[Vpe Upe]
% error('1')

%%%%%new translation
 

    del = K\-R(2:3*N+1);

    [Rratio max(abs(R(2:3*N+1)))]

    
    E(2:3:3*N-1) = eu(2:N+1,1);
    E(3:3:3*N) = eu(2:N+1,2);
    E(4:3:3*N+1) = eu(2:N+1,3);


     EE = E(2:3*N+1) + del;%*dt;
     E = NaN*ones(N+2,1);
     E(2:3*N+1) = EE;
     
     eu(2:N+1,1) = E(2:3:3*N-1) ;
     eu(2:N+1,2) = E(3:3:3*N) ;
     eu(2:N+1,3) = E(4:3:3*N+1);

%      if(count == 100)
%      eu
%      error('1')
%      end

%%%%%%new trans

    Upe = eu+u;


    for j = 2:N+1
%     e(j,1) = eu(j,1);
%     e(j,2) = eu(j,2);
%     e(j,3) = eu(j,3);
% [e(j,1),e(j,2),e(j,3)] = toprimitivevars(eu(j,1),eu(j,2),eu(j,3));
        [Vpe(j,1),Vpe(j,2),Vpe(j,3)] = toprimitivevars(Upe(j,1),Upe(j,2),Upe(j,3));
    end
     
    
    e = Vpe-V;
% error('1')
%%%%%%new trans


%      e
% exacterr

% error('1')
     
     t = t+dt;
     dtold = dt;
 end
 

    
    count

ee = e;
ee

w = exacterrv-ee

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
errerr2 = [sqrt(sum((exacterr(2:N+1,1)-ee(2:N+1,1)).^2)/N) sqrt(sum((exacterr(2:N+1,2)-ee(2:N+1,2)).^2)/N)  sqrt(sum((exacterr(2:N+1,3)-ee(2:N+1,3)).^2)/N) ]
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

