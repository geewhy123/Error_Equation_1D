function [uu,d] = irk1(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk2(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here


if (strcmp(eqn,'error')==1)
    val = NaN*ones(N+2,1);
    timesbet = 1*k+time;
    
    global TEND
    if(abs(timesbet(end)-TEND) < 1e-10)
        timesbet(end) = TEND;
    end
    
    for kk = 2:N+1
        val(kk,1:length(timesbet)) = fnval(timesbet,obj.Rsp(kk));
    end
    obj.errorSource = -1*val(:,1);
    
    
       istep = round(obj.curTime/obj.tStep) +1;
       -1*obj.residual(:,istep)
 obj.errorSource = -1*obj.residual(:,istep);
 
    
end

u0 = u;
d = 0;

unew = u;

Z = obj.unstructuredrecon(unew,p,eqn);



obj.curTime = obj.curTime + k;



f = 1;
uold = u;
while(max(abs(f)) > 1e-10)
    J = obj.computefluxjacobian((u),eqn);
    J = eye(N+2,N+2)-(k)*J;
    Z = obj.unstructuredrecon((u),p,eqn);
    f = u-uold-k*obj.computefluxintegral(Z,eqn);
    uold = u;
    u(2:N+1) = u(2:N+1)-J(2:N+1,2:N+1)\f(2:N+1);
end


uu = u;
uu(N+2) = NaN;

d = max(d,abs((uu-u0))/k);

end

