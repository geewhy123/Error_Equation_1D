function [uu,d] = rk1(u,x,f,k,h,N,p,t,phys);
%RK1 Summary of this function goes here
%   Detailed explanation goes here

uu = zeros(N+2,1);
[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
d = 0;
for i = 2:N+1

[upr,upl,delt] = reconflux(u,Z,f,k,h,i,N,p,phys);
%[delt]= updatesol(u,Z,f,k,h,i,N,p);
uu(i) = u(i)+k*delt;
d = max(d,abs(delt)); 
end
uu(N+2) = NaN;
   
end

