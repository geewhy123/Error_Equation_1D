function [ output_args ] = solvebyeulerjacobian( obj)
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

V = obj.initialSolution
U = NaN*ones(3*N+2,1);
obj.computeprimalpseudo();

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
%   if(obj.pOrder >= 4)
% % error('1')
%     dt = 0.00002;
%  end
%  kk = dt;
 Rold = R;
 dtold = 1;
 
 while(max(abs(R(2:3*N+1))) > 1e-11 )
     J = computeeulerfluxjacobian(obj,V,'solution');%,x,h,N,p);
    
     count = count +1;
         
     Rratio =norm(Rold(2:N+1),2)/norm(R(2:N+1),2); 
     dt = dtold*c2*Rratio;




 K = J(2:3*N+1,2:3*N+1)+eye(3*N)/dt;

[Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);

%  [er]=reconplot(x,h,N,p,Z);
Rold = R;
 [phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,'solution');%reconfluxsoln(Z,f,h,N,p,physics,t,obj)
 
R(2:3:3*N-1) = phi1(2:N+1);
R(3:3:3*N) = phi2(2:N+1);
R(4:3:3*N+1) = phi3(2:N+1);
u = NaN*ones(size(3*N+2,1));
for j = 2:N+1
[u(j,1),u(j,2),u(j,3)] = toconservedvars(V(j,1),V(j,2),V(j,3));
end




    del = K\-R(2:3*N+1);
    
%     if(mod(count,100)==0)
    max(abs(R(2:3*N+1)))
%     end

U(2:3:3*N-1) = u(2:N+1,1);
U(3:3:3*N) = u(2:N+1,2);
U(4:3:3*N+1) = u(2:N+1,3);


UU = U(2:3*N+1) + del;%*dt;
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

vv = u-ue;
%   cverr1 = sum(abs(vv(2:N+1)))/N
  cverr2 = sqrt(sum((vv(2:N+1)).^2)/N)

   plot(x,u,'*',x,ue,'o')


 obj.convSoln = u;
 
%  obj.convSoln

end

