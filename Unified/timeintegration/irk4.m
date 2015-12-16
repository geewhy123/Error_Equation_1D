function [uu,d] = irk4(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk2(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

c = [1/2+sqrt(3)/6 1/2-sqrt(3)/6];
if (strcmp(eqn,'error')==1)
    val = NaN*ones(N+2,1);
    timesbet = c*k+time;
    
    global TEND
    if(abs(timesbet(end)-TEND) < 1e-10)
        timesbet(end) = TEND;
    end
    
    for kk = 2:N+1
        val(kk,1:length(timesbet)) = fnval(timesbet,obj.Rsp(kk));
    end
%     obj.errorSource = -1*val(:,1:length(timesbet));
%     time
%     timesbet
%     val
%     error('1')
end

u0 = u;
d = 0;

unew = u;

% Z = obj.unstructuredrecon(unew,p,eqn);



% obj.curTime = obj.curTime + k;




uold = u;
% Z = obj.unstructuredrecon((u),p,eqn);
% y1 = obj.computefluxintegral(Z,eqn);
y1 =zeros(size(u));
y2 = y1;
a11 = 1/4;
a12 = 1/4+sqrt(3)/6;
a21 = 1/4-sqrt(3)/6;
a22 = a11;

f1 = zeros(N+2,1);
f2 = f1;
f = zeros(2*N+2,1);
f(2) = 1;
y = f;
J = zeros(2*N+2,2*N+2);
while(max(abs(f)) > 3e-9)
%     max(abs(f))


if (strcmp(eqn,'error')==1)
     obj.errorSource = -1*val(:,1);
     
         istep = round(obj.curTime/obj.tStep) +1;
         istep
 obj.errorSource = -1*obj.residual(:,istep);
end
    Z = obj.unstructuredrecon((uold+k*a11*y1+k*a12*y2),p,eqn);
    f1 = y1-obj.computefluxintegral(Z,eqn);
    J1 = obj.computefluxjacobian((uold+k*a11*y1+k*a12*y2),eqn);
    J11 = eye(N+2,N+2)-k*a11*J1;    
    J12 = -k*a12*J1;
  
    
%      q = 5*rand(size(f));
%      q1 = y1;
%      q2 = y2;
%      q1(2:N+1) = q(2:N+1);
%      q2(2:N+1) = q(N+2:2*N+1);
%      y1q = y1+q1;
%      y2q = y2+q2;
%      
%      Z = obj.unstructuredrecon((uold+k*a11*y1q+k*a12*y2q),p,eqn);
%     f1q = y1q-obj.computefluxintegral(Z,eqn);
% [f1q f1+J11*q1+J12*q2]
%     error('1')

if (strcmp(eqn,'error')==1)
     obj.errorSource = -1*val(:,2);
     
     
     
         istep = round(obj.curTime/obj.tStep) +1;
 obj.errorSource = -1*obj.residual(:,istep);
     
end


    Z = obj.unstructuredrecon((uold+k*a21*y1+k*a22*y2),p,eqn);
    f2 = y2-obj.computefluxintegral(Z,eqn);
    J2 = obj.computefluxjacobian((uold+k*a21*y1+k*a22*y2),eqn);
    J21 = -k*a21*J2;    
    J22 = eye(N+2,N+2)-k*a22*J2;
  
    
%      q1 = 5*rand(size(f1));
%      y1q = y1+q1;
%      Z = obj.unstructuredrecon((uold+k*a21*y1q+k*a22*y2q),p,eqn);
%     f2q = y2q-obj.computefluxintegral(Z,eqn);
%     [f2q f2+J21*q1+J22*q2]
%     error('1')
    
    
    J(2:2*N+1,2:2*N+1) = [J11(2:N+1,2:N+1) J12(2:N+1,2:N+1);J21(2:N+1,2:N+1) J22(2:N+1,2:N+1)];
    f(2:2*N+1) = [ f1(2:N+1);f2(2:N+1)];
        y(2:2*N+1) = [y1(2:N+1);y2(2:N+1)];
        
        
%         q = 0.5*rand(size(f));
%         y1q = y1;
%         y2q = y2;
%     y1q(2:N+1) = y1(2:N+1)+q(2:N+1);
%     y2q(2:N+1) = y2(2:N+1)+q(N+2:2*N+1);

%     Z = obj.unstructuredrecon((uold+k*a11*y1q+k*a12*y2q),p,eqn);
%     f1q = y1q-obj.computefluxintegral(Z,eqn);
%         Z = obj.unstructuredrecon((uold+k*a21*y1q+k*a22*y2q),p,eqn);
%     f2q = y2q-obj.computefluxintegral(Z,eqn);
  
% J
% J11
% J12
% J21
% J22
%         [[f1q(2:N+1);f2q(2:N+1)]-[[f(2:2*N+1)]+J(2:2*N+1,2:2*N+1)*q(2:2*N+1)]]
%         error('1')
        
%         size(f(2:2*N+1))
%         size(J(2:2*N+1,2:2*N+1))
%         size(y(2:2*N+1))
    y(2:2*N+1) = y(2:2*N+1)-J(2:2*N+1,2:2*N+1)\f(2:2*N+1);
     y1(2:N+1) = y(2:N+1);
    y2(2:N+1) = y(N+2:2*N+1);
    
    
  
%     error('1')
    max(abs(f));
end
 u(2:N+1) = u(2:N+1)+(k/2)*(y1(2:N+1)+y2(2:N+1));


uu = u;
uu(N+2) = NaN;

d = max(d,abs((uu-u0))/k);

end

