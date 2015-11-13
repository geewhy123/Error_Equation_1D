function [uu,d] = rk1(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)%rk1(u,x,f,k,h,N,p,t,phys);
%RK1 Summary of this function goes here
%   Detailed explanation goes here

uu = zeros(N+2,1);
% [Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
% 
Z = obj.unstructuredrecon(u,p,eqn);

d = 0;
% for i = 2:N+1
% 
% [upr,upl,delt] = reconflux(u,Z,f,k,h,i,N,p,phys);
% %[delt]= updatesol(u,Z,f,k,h,i,N,p);
% uu(i) = u(i)+k*delt;
% d = max(d,abs(delt)); 
% end
if(strcmp(eqn,'solution')==1)
    %     obj.curTime = obj.curTime + (steps-1)*k/length(c);
    phi = obj.computefluxintegral(Z,eqn);
elseif(strcmp(eqn,'error')==1)
    error('1')
%     obj.errorSource = -1*val(:,steps);
    phi(:,1) = obj.computefluxintegral(Z,'error');
else
    error('2')
end
uu = u+k*phi;

uu(N+2) = NaN;

 d = max(d,abs((uu-u))/k); 
end

