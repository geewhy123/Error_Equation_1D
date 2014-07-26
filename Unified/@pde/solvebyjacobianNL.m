function [errerr2,x,cverr2,exacterr,ee,te  ] = solvebyjacobianNL( obj )
%SOLVEBYJACOBIAN Summary of this function goes here
%   Detailed explanation goes here


ue = obj.exactSolution;
p = obj.pOrder;
q = obj.qOrder;
r = obj.rOrder;
h = obj.cellWidths;
N = obj.nCells;
x = obj.cellCentroids;
k = obj.tStep;
physics = obj.physics;
tlim = obj.endTime;

%%%%%%
% obj.hOrder = 0;
obj.computehigherpseudo();
%%%%%%

obj.computeprimalpseudo();






% 
% J = obj.computefluxjacobian(ue,'solution');%,x,h,N,p);
%  J
% error('1')

% [d,c] =  eig(J(2:N+1,2:N+1))
% nu = null(J(2:N+1,2:N+1));

%  error('1')
% J(2:N+1,2:N+1)
% 
% error('1')

[Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);



 % % % % % % % Z
% % % % % % % % % error('2')
% % % % % % % %   [er]=obj.reconplot(Z,'solution')%reconplot(x,h,N,p,Z)
% % % % % % % % %   error('1')
% % % % % % % %   
% % % % % % % %   
% % % % % % % % f = obj.source;
% % % % % % % %  [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)
% % % % % % % %  max(abs(tau))
% % % % % % % %  
% % % % % % % % %  sqrt(sum((tau(2:N+1)).^2)/N)
% % % % % % % % te = tau;
% % % % % % % % % 
% % % % % % % % % tau2 = tau;
% % % % % % % % % save('tau.mat','tau2','-append')
% % % % % % % % % error('1')
% % % % % % % % 
% % % % % % % % 
% % % % % % % % [sum(abs(tau(2:N+1)))/N sqrt(sum((tau(2:N+1).^2)/N)) max(abs(tau(2:N+1)))]
% % % % % % % % 
% % % % % % % % % errerr2= NaN;
% % % % % % % % % cverr2 = NaN;
% % % % % % % % % exacterr = NaN;
% % % % % % % % % ee = NaN;
% % % % % % % % % return;
% % % % % % % % figure
% % % % % % % % plot(x,tau,'o')
% % % % % % % % te=  sum(abs(tau(2:N+1)))/N;




obj.hOrder = obj.pOrder;
 refinecells = [];%[2 3 4 N-1 N N+1];%2 3 4 5 N-2 N-1 N N+1];
obj.refinecells = refinecells;
if(obj.hOrder > 0)
   obj.computehigherpseudo();
    [Zh] = obj.unstructuredrecon(ue,obj.hOrder,'solution'); 
    Znew = zeros(max(obj.hOrder,p),N+2);
    Znew(1:p,:) = Z;
    if(obj.hOrder >= p)
        for ii = 1:length(refinecells)
%     Znew(:,2) = Zh(:,2);
%     Znew(:,3) = Zh(:,3);
%     Znew(:,3) = Zh(:,4);
%     Znew(:,N) = Zh(:,N);
%     Znew(:,N+1) = Zh(:,N+1);
          Znew(:,refinecells(ii)) = Zh(:,refinecells(ii));
        end
    else
%     Znew(1:obj.hOrder,2) = Zh(1:obj.hOrder,2);
%     Znew(obj.hOrder+1:end,2) = zeros(p-obj.hOrder,1);
%     Znew(1:obj.hOrder,3) = Zh(1:obj.hOrder,3);
%     Znew(obj.hOrder+1:end,3) = zeros(p-obj.hOrder,1);
%     Znew(1:obj.hOrder,N) = Zh(1:obj.hOrder,N);
%     Znew(obj.hOrder+1:end,N) = zeros(p-obj.hOrder,1);
%     Znew(1:obj.hOrder,N+1) = Zh(1:obj.hOrder,N+1);
%     Znew(obj.hOrder+1:end,N+1) = zeros(p-obj.hOrder,1);
        for ii = 1:length(refinecells)
            Znew(1:obj.hOrder,refinecells(ii)) = Zh(1:obj.hOrder,refinecells(ii));
            Znew(obj.hOrder+1:end,refinecells(ii)) = zeros(p-obj.hOrder,1);
        end
    end
    Znew
%     error('1')
    [tau]=obj.computefluxintegral(Znew,'solution');
end

[sum(abs(tau(2:N+1)))/N sqrt(sum((tau(2:N+1).^2)/N)) max(abs(tau(2:N+1)))]

figure
plot(x,tau,'o')

figure
% obj.pOrder = obj.hOrder;
obj.reconplot(Znew,'solution')
% obj.reconplot(Z,'solution')
% error('1') 
Znew
te = tau;






%   [er]=obj.reconplot(Z,'solution')%x,h,N,p,Z);
%   error('1')
f = obj.source;
%  [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

%   tau
te1 = sum(abs(tau(2:N+1)))/N 

te = tau;


% tau6 = tau;
% save('tauB.mat','tau6','-append')
% error('1')


% errerr2= NaN;
% cverr2 = NaN;
% exacterr = NaN;
% ee = NaN;
% return;


%   error('1')
 
 max(abs(tau))
  
 
 
J = obj.computefluxjacobian(ue,'solution');%,x,h,N,p);
 
 
 del = ones(N,1);
 R = ones(N+2,1);
 t=0;
 u = ue;
 count = 0;
 c2 = 10;
 dt = 0.001;
%   if(obj.pOrder >= 4)
% % error('1')
%     dt = 0.00002;
%  end
%  kk = dt;
 Rold = R;
 dtold = 1;
 while(max(abs(R)) > 1e-11 )
     J = obj.computefluxjacobian(u,'solution');%,x,h,N,p);
    
     count = count +1;
         
     Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 
     dt = dtold*c2*Rratio;




 K = J(2:N+1,2:N+1)+eye(N)/dt;
% K
% error('1')
%  K = (K+K')/2;

[Z] = obj.unstructuredrecon(u,p,'solution');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
Rold = R;
 [R]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
    del = K\-R(2:N+1);
    
%     if(mod(count,100)==0)
    max(abs(R(2:N+1)))
%     end


     uu = u(2:N+1) + del;%*dt;
     u = NaN*ones(N+2,1);
     u(2:N+1) = uu;
     t = t+dt;
     
     
     dtold = dt;
     
 end

 
%   u
%   error('1')
%   max(abs(u-ue))
vv = ue-u;
%   cverr1 = sum(abs(vv(2:N+1)))/N
  cverr2 = sqrt(sum((vv(2:N+1)).^2)/N)

   plot(x,u,'*',x,ue,'o')
%    error('1')
%   u-ue

%  v=J(2:N+1,2:N+1)\tau(2:N+1)
% 
%   if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(J(2:N+1,2:N+1)))) < 1e-5)
%     v = pinv(J(2:N+1,2:N+1))*tau(2:N+1);
%   end
% %  max(abs(v-ue(2:N+1)))
% max(abs(v))
% %  error('1')
%  figure
%  plot(x,u-ue,x(2:N+1),v)
%  
%  v(2:N+1) = v;
%  v(1) = NaN;
%  v(N+2) = NaN;
%  
%  u = ue-v
%  plot(x,u)
%  v
% %  ue-v
% % cverr1 = sum(abs(v(2:N+1)))/N
% %  cverr2 = sqrt(sum((v(2:N+1)).^2)/N)

 obj.convSoln = u;
 u
vv
%  obj.convSoln
%  error('1')
 count
% error('1')





% % % % accuracy
% % % % pp=obj.pOrder;
% % % % obj.pOrder = obj.qOrder;
% % % % p = obj.pOrder;
% % % % s = f;
% % % % obj.hOrder = 0;obj.pOrder;
% % % % obj.computeprimalpseudo();
% % % %  [Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  [tauq]=obj.computefluxintegral(Z,'solution');
% % % % 
% % % %  tauq
% % % %  obj.source = obj.source*0;
% % % %   [Z] = obj.unstructuredrecon(u,q,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  FIq = obj.computefluxintegral(Z,'solution');
% % % %  FIq
% % % % %  obj.computerespseudo();
% % % % obj.pOrder = r;
% % % %  obj.computeprimalpseudo();
% % % %   [Z] = obj.unstructuredrecon(u,r,'solution');%ue,x,h,N,NaN,NaN,p);
% % % %  FIr = obj.computefluxintegral(Z,'solution'); 
% % % %   FIr
% % % %  
% % % % 
% % % % % mean(abs(tauq(2:N+1)))
% % % % g=tauq-(FIq-FIr)
% % % % 
% % % % J = obj.computefluxjacobian(ue,'solution');%,x,h,N,p);
% % % % J(2:N+1,2:N+1)\(tauq(2:N+1)-FIq(2:N+1)+s(2:N+1))
% % % % J(2:N+1,2:N+1)\(-FIr(2:N+1)+s(2:N+1))
% % % % ue-u
% % % % % -FIr+s
% % % % mean(abs(g(2:N+1)))
% % % % 
% % % % b = tauq+s-FIq
% % % % mean(abs(b(2:N+1)))
% % % % c = FIr-s
% % % % mean(abs(c(2:N+1)))
% % % % 
% % % %  error('1')
% % % % obj.pOrder = pp;
% % % % p = pp;
% % % % obj.source = s;
% % % % 
 



%%%%residual

if(q> 0 && r>0)

  obj.computerespseudo();
  
  [Zr] = obj.unstructuredrecon(u,r,'residual');
  figure
  obj.reconplot(Zr,'residual');
  
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
Rend = obj.computefluxintegral(Zr,'residual');
 
Rend
R1  = sum(abs(Rend(2:N+1)))/N
%    error('2')
 
 

% clearvars -except x obj 
obj.computeerrorpseudo();


 Zu = obj.unstructuredrecon(obj.convSoln,obj.qOrder,'error');
 obj.convSolnRecon = Zu;
 

% test using TE for p~=q 
% load('tauB.mat')
% [Zq] = obj.unstructuredrecon(u,q,'error');%ue,x,h,N,NaN,NaN,p);
% Rq = obj.computefluxintegral(Zq,'residual');

% [tau4 tau4-Rq -Rend]
% error('1')



%%%%error equation

if(obj.bcLeftType == 'D')
   obj.bcLeftVal = 0; 
end
if(obj.bcRightType == 'D')
    obj.bcRightVal = 0;
end

exacterr = ue-u;

%  obj.computeerrorpseudo();
[Z] = obj.unstructuredrecon(exacterr,q,'error');%ue,x,h,N,NaN,NaN,p);

% Z
% figure
% obj.reconplot(Z,'error')
% hold on
% plot(x,exacterr)
% error('1')


f = -Rend;
   obj.errorSource = f;%tau6-Rq;%f;%tau;
 
%    [f -FIr+s]
%    error('1')
 
 Je = obj.computefluxjacobian(exacterr,'error');
% obj.errorRM
% error('1')



% error('1')
%  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,obj)
%  error('1')
 [tauE]= obj.computefluxintegral(Z,'error')
%   error('1')
 
% error('1')

 del = ones(N,1);
 R = ones(N+2,1);
 t=0;
 
 e = exacterr
 ee = ones(N+2,1);
 ee(1)=NaN;
 ee(N+2) = NaN;
 count = 0;
 
  dt = 0.001;
%  if(obj.qOrder > 4)
% % error('1')
%     kk = 0.00005;
%  end
  
  
 while(max(abs(R)) > 1e-11 )
     Je = obj.computefluxjacobian(e,'error');%,x,h,N,p);
    
     
%      Jue = obj.computefluxjacobian(u+e,'error');%,x,h,N,p);
%      Ju = obj.computefluxjacobian(u,'error');%,x,h,N,p);
%      Je = Jue-Ju;
     
     
     
     count = count +1;
%      if(count < 50)
%         dt = kk*(40/N)^2; 
        
      Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 

           dt = dtold*c2*Rratio;



 K = Je(2:N+1,2:N+1)+eye(N)/dt;


[Z] = obj.unstructuredrecon(e,q,'error');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
Rold = R;
 [R]=obj.computefluxintegral(Z,'error');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
 
 
   [K\-R(2:N+1) pinv(K)*-R(2:N+1) (K'*K)\(K'*-R(2:N+1)) R(2:N+1)]
%   error('1')
    del = (K'*K)\(K'*-R(2:N+1));%pinv(K)*-R(2:N+1);%K\-R(2:N+1);
    
%     if(mod(count,100)==0)
    max(abs(R(2:N+1)))
%     end
    
     ee(2:N+1) = e(2:N+1) + del;%*dt;
     e = NaN*ones(N+2,1);
     e = ee;
     t = t+dt;
     dtold = dt;
 end



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
plot(x,ee,'*',x,exacterr,'o')

-tauE
ee
else
   errerr2 = NaN;
   exacterr = NaN;
   ee = NaN;
    
end


% J(2:N+1,2:N+1)*exacterr(2:N+1)
% tau(2:N+1)
%  J(2:N+1,2:N+1)*ue(2:N+1)
%  obj.source(2:N+1)
% Je(2:N+1,2:N+1)*w(2:N+1)
% tauE(2:N+1)
% error('1')
end

