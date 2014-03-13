function [uu,d] = rk7(u,x,f,k,h,N,p,t,phys,uder,time,gsp)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

% b = [7/90 0 32/90 12/90 32/90 7/90];
% A=[

%c = [1/6 1/3 1/2 2/11 2/3 6/7 0 1]

phi = zeros(N+2,1);
phiII = zeros(N+2,1);
phiIII = zeros(N+2,1);
phiIV = zeros(N+2,1);
phiV = zeros(N+2,1);
phiVI = zeros(N+2,1);
phiVII = zeros(N+2,1);
phiVIII = zeros(N+2,1);
phiIX = zeros(N+2,1);
uII = zeros(N+2,1);
uIII = zeros(N+2,1);
uIV = zeros(N+2,1);
uV = zeros(N+2,1);
uVI = zeros(N+2,1);
uVII = zeros(N+2,1);
uVIII = zeros(N+2,1);
uIX = zeros(N+2,1);
uu = zeros(N+2,1);

[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
d = 0;
for i = 2:N+1
       
[upr,upl,phi(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time,gsp);
uII(i) = u(i)+(k/6)*phi(i);
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k/6,gsp);
uIII(i) = u(i)+(k/3)*(phiII(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uIII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k/3,gsp);
uIV(i) = u(i)+(k/8)*(phi(i)+3*phiIII(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uIV,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIV(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k/2,gsp);
uV(i) = u(i)+(k/1331)*(148*phi(i)+150*phiIII(i)-56*phiIV(i));
   
end

[Z]=unstructuredrecon(uV,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiV(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k*(2/11),gsp);
uVI(i) = u(i)+(k/1701)*(-2828*phi(i)-10710*phiIII(i)+4024*phiIV(i)+10648*phiV(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uVI,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiVI(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k*(2/3),gsp);
uVII(i) = u(i)+(k/16807)*(17262*phi(i)+60858*phiIII(i)-19176*phiIV(i)-51909*phiV(i)+7371*phiVI(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uVII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiVII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k*(6/7),gsp);
uVIII(i) = u(i)+(k)*((5/154)*phi(i)+(96/539)*phiIV(i)-(1815/20384)*phiV(i)-(405/2464)*phiVI(i)+(49/1144)*phiVII(i));
   
%d = max(d,abs(delt)); 
end

[Z]=unstructuredrecon(uVIII,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiVIII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time,gsp);
uIX(i) = u(i)+(k)*((-113/32)*phi(i)-(195/22)*phiIII(i)+(32/7)*phiIV(i)+(29403/3584)*phiV(i)-(729/512)*phiVI(i)+(1029/1408)*phiVII(i)+(21/16)*phiVIII(i));
   
%d = max(d,abs(delt)); 
end


[Z]=unstructuredrecon(uIX,x,h,N,NaN,NaN,p);
for i = 2:N+1        
       
[upr,upl,phiIX(i)] = reconflux(u,Z,f,k,h,i,N,p,phys,uder,NaN,time+k,gsp);

uu(i) = u(i)+(k)*((32/105)*phiIV(i)+(1771561/6289920)*phiV(i)+(243/2560)*phiVI(i)+(16807/74880)*phiVII(i)+(77/1440)*phiVIII(i)+(11/270)*phiIX(i));   
d = max(d,abs(((32/105)*phiIV(i)+(1771561/6289920)*phiV(i)+(243/2560)*phiVI(i)+(16807/74880)*phiVII(i)+(77/1440)*phiVIII(i)+(11/270)*phiIX(i)))); 
end
uu(N+2) = NaN;
   
end

