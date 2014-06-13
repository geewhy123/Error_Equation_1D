function [ output_args ] = computeeulerjacobian( input_args )
%COMPUTEEULERJACOBIAN Summary of this function goes here
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

[Z] = obj.unstructuredrecon(ue,p,'solution');%ue,x,h,N,NaN,NaN,p);

%   [er]=obj.reconplot(Z,'solution')%x,h,N,p,Z);
%   error('1')
f = obj.source;
 [tau]=obj.computefluxintegral(Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,tlim,obj)

  tau
te1 = sum(abs(tau(2:N+1)))/N 

te = tau;


% errerr2= NaN;
% cverr2 = NaN;
% exacterr = NaN;
% ee = NaN;
% return;


%   error('1')
 
 max(abs(tau))
  
 del = ones(N,1);
 R = ones(N+2,1);
 t=0;
 u = ue;
 count = 0;
 c2 = 10;
 dt = 0.0001;
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
         
%      Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 
%      dt = dtold*c2*Rratio;




 K = J(2:N+1,2:N+1)+eye(N)/dt;

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
vv = u-ue;
%   cverr1 = sum(abs(vv(2:N+1)))/N
  cverr2 = sqrt(sum((vv(2:N+1)).^2)/N)

   plot(x,u,'*',x,ue,'o')


 obj.convSoln = u;
 
%  obj.convSoln

end

