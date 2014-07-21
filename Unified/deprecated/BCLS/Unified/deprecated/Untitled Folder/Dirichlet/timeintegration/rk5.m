function [uu,d] = rk5(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

% b = [7/90 0 32/90 12/90 32/90 7/90];
% A=[


uu = zeros(N+2,1);
[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
d = 0;
for i = 2:N+1
       
[upr,upl,phi(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uII(i) = u(i)+(k/4)*phi(i);
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uIII(i) = u(i)+(k/8)*(phi(i)+phiII(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uIII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uIV(i) = u(i)+(k/2)*phiIII(i);
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uIV,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIV(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uV(i) = u(i)+(k/16)*(3*phi(i)-6*phiII(i)+6*phiIII(i)+9*phiIV(i));
   
end

[Z]=unstructuredrecon(uV,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiV(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
uVI(i) = u(i)+(k/7)*(-3*phi(i)+8*phiII(i)+6*phiIII(i)-12*phiIV(i)+8*phiV(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uVI,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiVI(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);

uu(i) = u(i)+(k/90)*((7)*phi(i)+(0)*phiII(i)+(32)*phiIII(i)+(12)*phiIV(i)+(32)*phiV(i)+(7)*phiVI(i));   
d = max(d,abs(((7)*phi(i)+(0)*phiII(i)+(32)*phiIII(i)+(12)*phiIV(i)+(32)*phiV(i)+(7)*phiVI(i))/90)); 
end
uu(N+2) = NaN;
   
end

