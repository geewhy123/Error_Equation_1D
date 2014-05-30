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
% J
% error('1')

% [d,c] =  eig(J(2:N+1,2:N+1))


% nu = null(J(2:N+1,2:N+1));

%  error('1')
% J(2:N+1,2:N+1)
% 
% error('1')

[Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
f = obj.source;
 [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

 tau
%  error('1')
 
 max(abs(tau))
  
 del = ones(N,1);
 t=0;
 u0 = ue;
 count = 0;
 while(max(abs(del)) > 1e-10 )
     count = count +1;
     dt = .01;
     
 K = J(2:N+1,2:N+1)+eye(N)/dt;
% K
% error('1')
%  K = (K+K')/2;

[Z] = obj.unstructuredrecon(u0,p,'solution');%u0,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
 [R]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
    del = K\R(2:N+1);
    max(abs(R(2:N+1)))
    
     uu = u0(2:N+1) + del;%*dt;
     u0 = NaN*ones(N+2,1);
     u0(2:N+1) = uu;
     t = t+dt;
 end
 
 
  u0
  error('1')
%  max(abs(u0-ue))
%  u0-ue

 v=J(2:N+1,2:N+1)\tau(2:N+1)

  if(obj.bcLeftType=='P' && obj.bcRightType == 'P' && min(abs(eig(J(2:N+1,2:N+1)))) < 1e-5)
    v = pinv(J(2:N+1,2:N+1))*tau(2:N+1);
  end
%  max(abs(v-ue(2:N+1)))
max(abs(v))
%  error('1')
 figure
 plot(x,u0-ue,x(2:N+1),v)
 
 v(2:N+1) = v;
 v(1) = NaN;
 v(N+2) = NaN;
 
 u = ue-v
 plot(x,u)
 v
%  ue-v
 cverr2 = sqrt(sum((v(2:N+1)).^2)/N)

 
 
error('1')




  obj.computerespseudo();
  [Zr] = obj.unstructuredrecon(u,r,'residual');
%     [left,right] = computeflux(Zr,h,N,r,physics,'residual',obj);
%     Rend= (right-left)./h-f;
Rend = obj.computefluxintegral(Zr,'residual');
 
 
 
 
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

