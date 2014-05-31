function [errerr2,x,cverr2,exacterr,ee  ] = solvebyjacobianNL( obj )
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


obj.computeprimalpseudo();

J = obj.computefluxjacobian(ue,'solution');%,x,h,N,p);
%  J
% error('1')

% [d,c] =  eig(J(2:N+1,2:N+1))
% nu = null(J(2:N+1,2:N+1));

%  error('1')
% J(2:N+1,2:N+1)
% 
% error('1')

[Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);

%  [er]=obj.reconplot(Z)%x,h,N,p,Z);
%  error('1')
f = obj.source;
 [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

 tau
%  error('1')
 
 max(abs(tau))
  
 del = ones(N,1);
 R = ones(N+2,1);
 t=0;
 u = ue;
 count = 0;
 while(max(abs(R)) > 1e-11 )
%      J = obj.computefluxjacobian(u,'solution');%,x,h,N,p);
    
     count = count +1;
%      if(count < 50)
        dt = 0.001; 
        
%      elseif(count <100)
%          dt = 0.0001;
%      elseif(count < 500)
%          dt = 0.0005;
%      elseif(count < 1000)
%          dt = 0.001;
%      elseif(count < 10000)
%          dt = 0.01;
%      else
%          dt = 0.05;
%      end
     
% if (mod(count,100)==0)
%    dt = dt*(count/100);
%    error('3')
% end


 K = J(2:N+1,2:N+1)+eye(N)/dt;
% K
% error('1')
%  K = (K+K')/2;

[Z] = obj.unstructuredrecon(u,p,'solution');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
 [R]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
    del = K\R(2:N+1);
    max(abs(R(2:N+1)))
    
     uu = u(2:N+1) + del;%*dt;
     u = NaN*ones(N+2,1);
     u(2:N+1) = uu;
     t = t+dt;
 end
 
 
%   u
%   error('1')
%   max(abs(u-ue))
vv = u-ue;
%   cverr1 = sum(abs(vv(2:N+1)))/N
  cverr2 = sqrt(sum((vv(2:N+1)).^2)/N)

   plot(x,u,'*')
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

 
 count
% error('1')




  obj.computerespseudo();
  [Zr] = obj.unstructuredrecon(u,r,'residual');
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
Rend = obj.computefluxintegral(Zr,'residual');
 
Rend
R2 = sqrt(sum((Rend(2:N+1)).^2)/N)
 error('2')
 
 
 obj.computeerrorpseudo();
[Z] = obj.unstructuredrecon(ue-u,q,'error');%ue,x,h,N,NaN,NaN,p);
f = -Rend;
   obj.errorSource = f;
%    f
%    error('2')
 
 
 Je = obj.computefluxjacobian(ue,'error');
% obj.errorRM
% error('1')



% error('1')
%  [tauE]=reconfluxsoln(Z,f,h,N,q,physics,tlim,obj)
%  error('1')
 [tauE]= obj.computefluxintegral(Z,'error')
%  error('1')
 
%  Je
%  error('1')
w = Je(2:N+1,2:N+1)\tauE(2:N+1)

% x1 = ones(N,1);
%  null(Je(2:N+1,2:N+1),'r')
% Je
%   error('1')
if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(Je(2:N+1,2:N+1)))) < 1e-5)
%   error('5')
%   [Q,R] = qr(Je(2:N+1,2:N+1)') ;
% w = Q*(R'\tauE(2:N+1)) ;
% w = Je(2:N+1,2:N+1)'*Je(2:N+1,2:N+1)\Je(2:N+1,2:N+1)'*tauE(2:N+1)

% w = w-dot(w,x1')*x1*dot(w,w);
    w = pinv(Je(2:N+1,2:N+1))*tauE(2:N+1)
  end


max(abs(w))
% error('1')
figure
plot(x(2:N+1),w)

% error('2')
 

w(2:N+1) = w;
w(1) = NaN;
w(N+2) = NaN;
% cverr2 = sqrt(sum((ue(2:N+1)-u(2:N+1)).^2)/N);
exacterr = ue-u;
ee = exacterr - w;
errerr2 = sqrt(sum((exacterr(2:N+1)-ee(2:N+1)).^2)/N) 
w
 
end

