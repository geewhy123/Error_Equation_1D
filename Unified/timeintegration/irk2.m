function [uu,d] = irk2(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk2(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here


if (strcmp(eqn,'error')==1)
    val = NaN*ones(N+2,1);
    timesbet = 0.5*k+time;
    
    global TEND
    if(abs(timesbet(end)-TEND) < 1e-10)
%         error('1')
        timesbet(end) = TEND;
    end
    
    for kk = 2:N+1
        val(kk,1:length(timesbet)) = fnval(timesbet,obj.Rsp(kk));
        uval(kk,1:length(timesbet)) = fnval(timesbet,obj.Usp(kk)); 
        utval(kk,1:length(timesbet)) = fnval(fnder(obj.Usp(kk)),timesbet); 
        
        
        
    end
    obj.errorSource = -1*val(:,1);
    
%      obj.computerespseudo();
    Z = obj.unstructuredrecon(uval,obj.rOrder,'residual');
    r = obj.computefluxintegral(Z,'residual');
utval(N+2) = NaN;

% %     obj.errorSource = (utval-r);

% 
% i = round(obj.curTime/obj.tStep) +1
%     obj.errorSource = obj.Rall(:,2*i);
%     [obj.Rall(:,2*i) -1*val(:,1)]
% %     error('1')

% [(utval-r) -1*val(:,1)]
  
% error('1')
    
end

u0 = u;
d = 0;

unew = u;

Z = obj.unstructuredrecon(unew,p,eqn);



% obj.curTime = obj.curTime + k;


% Z = obj.unstructuredrecon(u,p,eqn);
% f = obj.computefluxintegral(Z,eqn);
% J = obj.computefluxjacobian(u,eqn);
% [u f]
% d= 0.00025*rand(size(u));
% up  = u+d;
% Z = obj.unstructuredrecon(up,p,eqn);
% fp = obj.computefluxintegral(Z,eqn);
% [up fp]
% [fp-f-J*d]
% max(abs(fp-f-J*d))
% error('1')
obj.curTime;
f = 1;
uold = u;
y = u;
while(max(abs(f)) > 1e-9)
%     J = obj.computefluxjacobian((uold+u)/2,eqn);
%     J = eye(N+2,N+2)-(k/2)*J;
%     Z = obj.unstructuredrecon((uold+u)/2,p,eqn);
%     f = u-uold-k*obj.computefluxintegral(Z,eqn);
%         uold = u;
%     u(2:N+1) = u(2:N+1)-J(2:N+1,2:N+1)\f(2:N+1);
    

    
    
    Z = obj.unstructuredrecon(u+k*y/2,p,eqn);
    f = y-obj.computefluxintegral(Z,eqn);
    J = obj.computefluxjacobian(u+k*y/2,eqn);
    J = eye(N+2,N+2)-(k/2)*J;
    y(2:N+1) = y(2:N+1)-J(2:N+1,2:N+1)\f(2:N+1);
    max(abs(f));    
end

u(2:N+1) = u(2:N+1)+k*y(2:N+1);
% obj.curTime = obj.curTime + k;

uu = u;
uu(N+2) = NaN;

d = max(d,abs((uu-u0))/k);


if(strcmp(eqn,'solution')==1)

i = round(obj.curTime/obj.tStep) +1;

    obj.Rall(:,2*i) = y;
     Z = obj.unstructuredrecon(u,p,eqn);
     f = obj.computefluxintegral(Z,eqn);
    obj.Rall(:,2*i+1) = f;
end
end

