function [uu,d] = rk4(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

uu = zeros(N+2,1);
[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
d = 0;
for i = 2:N+1
       
[upr,upl,phi(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uI(i) = u(i)+(k/2)*phi(i);
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uI,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiI(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uII(i) = u(i)+(k/2)*phiI(i);
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uIII(i) = u(i)+(k)*phiII(i);
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uIII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);

uu(i) = u(i)+(k/6)*(phi(i)+2*phiI(i)+2*phiII(i)+phiIII(i));   
d = max(d,abs((phi(i)+2*phiI(i)+2*phiII(i)+phiIII(i))/6)); 
end
uu(N+2) = NaN;
   
end

