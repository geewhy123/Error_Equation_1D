function [uu,d] = rk2(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk2(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here
c = [0 1/2 ];


if (strcmp(eqn,'error')==1)
    val = NaN*ones(N+2,length(c));
    timesbet = c*k+time;
    
    global TEND
    if(abs(timesbet(end)-TEND) < 1e-10)
        timesbet(end) = TEND;
    end
    
    for kk = 2:N+1
        val(kk,1:length(timesbet)) = fnval(timesbet,obj.Rsp(kk));
    end
end
uu = zeros(N+2,1);
Z = obj.unstructuredrecon(u,p,eqn);

d = 0;

A = [1/2 0; 0 1];
    unew = u;
     phi = zeros(N+2,length(c));
for steps = 1:length(c)
Z = obj.unstructuredrecon(unew,p,eqn);
    if(strcmp(eqn,'solution')==1)
            obj.curTime = obj.curTime + (steps-1)*k/length(c);
            
        phi(:,steps) = obj.computefluxintegral(Z,eqn);
      
     

    elseif(strcmp(eqn,'error')==1)
        obj.errorSource = -1*val(:,steps);
        %         obj.curTime = obj.curTime + (steps-1)*k/length(c);
        phi(:,steps) = obj.computefluxintegral(Z,'error');

    else
        error('2')
    end
    unew = u+phi*A(steps,:)'*k;
   
   
end
% phi
% error('1')
% unew
% error('2')
    uu = unew;
uu(N+2) = NaN;

 d = max(d,abs((uu-u))/k); 

end

