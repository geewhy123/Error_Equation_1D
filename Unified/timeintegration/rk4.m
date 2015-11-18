function [uu,d] = rk4(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk4(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

% uu = zeros(N+2,1);
% [Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
% d = 0;
% for i = 2:N+1
%        
% [upr,upl,phi(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
% uI(i) = u(i)+(k/2)*phi(i);
% %d = max(d,abs(delt)); 
% end
% 
% [Z]=unstructuredrecon(uI,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
%        
% [upr,upl,phiI(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
% uII(i) = u(i)+(k/2)*phiI(i);
%    
% %d = max(d,abs(delt)); 
% end
% 
% [Z]=unstructuredrecon(uII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
%        
% [upr,upl,phiII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
% uIII(i) = u(i)+(k)*phiII(i);
%    
% %d = max(d,abs(delt)); 
% end
% 
% [Z]=unstructuredrecon(uIII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
%        
% [upr,upl,phiIII(i)] = reconflux(u,Z,f,k,h,i,N,p,phys);
% 
% uu(i) = u(i)+(k/6)*(phi(i)+2*phiI(i)+2*phiII(i)+phiIII(i));   
% d = max(d,abs((phi(i)+2*phiI(i)+2*phiII(i)+phiIII(i))/6)); 
% end
% uu(N+2) = NaN;
%    

uu = zeros(N+2,1);
Z = obj.unstructuredrecon(u,p,eqn);

d = 0;
c = [0 1/2 1/2 1];
A = [1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 1/3 1/3 1/6];
    unew = u;
     phi = zeros(N+2,length(c));
for steps = 1:length(c)
Z = obj.unstructuredrecon(unew,p,eqn);
    if(strcmp(eqn,'solution')==1)
            obj.curTime = obj.curTime + (steps-1)*k/length(c);
            
        phi(:,steps) = obj.computefluxintegral(Z,eqn);
     

    elseif(strcmp(eqn,'error')==1)
        error('1')
        %     obj.errorSource = -1*val(:,steps);
        phi(:,1) = obj.computefluxintegral(Z,'error');
    else
        error('2')
    end
    unew = u+phi*A(steps,:)'*k;
   
   
end

    uu = unew;
uu(N+2) = NaN;

 d = max(d,abs((uu-u))/k); 


end

