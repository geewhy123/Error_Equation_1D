function [uu,d] = icn(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)% rk2(u,x,f,k,h,N,p,t,phys)
%RK1 Summary of this function goes here
%   Detailed explanation goes here


if (strcmp(eqn,'error')==1)
    val0 = NaN*ones(N+2,1);
    val1 = NaN*ones(N+2,1);
    
    for kk = 2:N+1
        val0(kk,1) = 999*fnval(time,obj.Rsp(kk));
        val1(kk,1) = 999*fnval(time+k,obj.Rsp(kk));
    end
    obj.errorSource = 0*val0(:,1);
    i = round(obj.curTime/obj.tStep) +1;
    R0 = obj.residual(:,i);
    R1 = obj.residual(:,i+1);
%     
%        J = obj.primalJacobian;

    J0 = obj.computefluxjacobian(obj.Uall(:,i),'solution');
    JJ0 = J0(2:N+1,2:N+1);

    J1 = obj.computefluxjacobian(obj.Uall(:,i+1),'solution');
    JJ1 = J1(2:N+1,2:N+1);
    
    
    Z0 = obj.unstructuredrecon(obj.Uall(:,i),obj.pOrder,'solution');
    f0 = obj.computefluxintegral(Z0,'solution');
    Z1 = obj.unstructuredrecon(obj.Uall(:,i+1),obj.pOrder,'solution');
    f1 = obj.computefluxintegral(Z1,'solution');
     
% %     R0(2:N+1) =  -(1/12)*obj.tStep^2*JJ^3*obj.Uall(2:N+1,i);%+(1/12)*obj.tStep^3*JJ^4*obj.Uall(2:N+1,i);
% %     R1(2:N+1) =  -(1/12)*obj.tStep^2*JJ^3*obj.Uall(2:N+1,i+1);%+(1/12)*obj.tStep^3*JJ^4*obj.Uall(2:N+1,i+1);
% 
%      R0(2:N+1) =  -(1/12)*obj.tStep^2*JJ0^2*f0(2:N+1);
%     R1(2:N+1) =  -(1/12)*obj.tStep^2*JJ1^2*f1(2:N+1);
%     
% %     R0(2:N+1) = R0(2:N+1) + -(1/24)*obj.tStep^3*JJ0^3*f0(2:N+1);
% %     R1(2:N+1) = R1(2:N+1) +  -(1/24)*obj.tStep^3*JJ1^3*f1(2:N+1);
% 
% H0 = zeros(N+2,N+2,N+2);
% H1 = H0;
% for kk = 2:N+1
%     delta = 1e-8;
%     del = zeros(N+2,1);
%     del(kk) = delta;
%    H0(:,:,kk) =  (obj.computefluxjacobian(obj.Uall(:,i)+del,'solution')-obj.computefluxjacobian(obj.Uall(:,i),'solution'))/delta;
%    H1(:,:,kk) =  (obj.computefluxjacobian(obj.Uall(:,i+1)+del,'solution')-obj.computefluxjacobian(obj.Uall(:,i+1),'solution'))/delta;
% end
% for ii = 2:N+1
%     for jj = 2:N+1
%         for kk = 2:N+1
%             R0(ii) = R0(ii)+H0(ii,jj,kk)*f0(jj)*f0(kk)*(-1/12)*obj.tStep^2;
%             R1(ii) = R1(ii)+H1(ii,jj,kk)*f1(jj)*f1(kk)*(-1/12)*obj.tStep^2;
%         end
%     end 
% end


% % grad(Jf)f
% R0 = R0*0;
% R1 = R1*0;
%     delta = 1e-8;
% for kk = 2:N+1
%     Zd = obj.unstructuredrecon(obj.Uall(:,i)+delta*f0,obj.pOrder,'solution');
%     fd = obj.computefluxintegral(Zd,'solution');
%     Jd = obj.computefluxjacobian(obj.Uall(:,i)+delta*f0,'solution');
%     R0(2:N+1) =  ((Jd(2:N+1,2:N+1)*fd(2:N+1)-J0(2:N+1,2:N+1)*f0(2:N+1))/delta)*(-1/12)*obj.tStep^2;
%     
%     Zd = obj.unstructuredrecon(obj.Uall(:,i+1)+delta*f1,obj.pOrder,'solution');
%     fd = obj.computefluxintegral(Zd,'solution');
%     Jd = obj.computefluxjacobian(obj.Uall(:,i+1)+delta*f1,'solution');
%     R1(2:N+1) =  ((Jd(2:N+1,2:N+1)*fd(2:N+1)-J1(2:N+1,2:N+1)*f1(2:N+1))/delta)*(-1/12)*obj.tStep^2;
% 
% end

    
% unsteady burgers?
% if(i==4)
%    plot(R0+R1)
%    error('1')
% end
    
%     ZZ0 = obj.unstructuredrecon(obj.Uall(:,i),p,'error');
%     ff0 = obj.computefluxintegral(ZZ0,'solution');
%     R0 = (obj.Uall(:,i+1)-obj.Uall(:,i))/obj.tStep-ff0;
%      ZZ1 = obj.unstructuredrecon(obj.Uall(:,i+1),p,'error');
%     ff1 = obj.computefluxintegral(ZZ1,'solution');
%     R1 = (obj.Uall(:,i+1)-obj.Uall(:,i))/obj.tStep-ff1;
    
    
    
%   [R0 R1]
% figure
%  plot(x,R0,'*')
%  error('1')
 


% %%%
% R0 = 0*R0;
% R1 = 0*R1;
% Z0 = obj.unstructuredrecon(obj.Uall(:,i),obj.pOrder,'solution');
% y1 = k*obj.computefluxintegral(Z0,'solution');
% Z0 = obj.unstructuredrecon(obj.Uall(:,i)+y1/2,obj.pOrder,'solution');
% y2 = k*obj.computefluxintegral(Z0,'solution');
% Z0 = obj.unstructuredrecon(obj.Uall(:,i)+y2/2,obj.pOrder,'solution');
% y3 = k*obj.computefluxintegral(Z0,'solution');
% Z0 = obj.unstructuredrecon(obj.Uall(:,i)+y3,obj.pOrder,'solution');
% y4 = k*obj.computefluxintegral(Z0,'solution');
% R0 = (1/1)*(obj.Uall(:,i+1)-obj.Uall(:,i))-(1/1)*(y1/6-y2/3-y3/3-y4/6);
% R1 = R0;
%%%

%     %exact te
%     r0 = R0;r1 = R1;
%     if(obj.pOrder == obj.qOrder)
%     ZZ1 = obj.unstructuredrecon(obj.exactSolutionAll(:,i+1),p,'solution');
%     ff1 = obj.computefluxintegral(ZZ1,'solution');
%     ZZ0 = obj.unstructuredrecon(obj.exactSolutionAll(:,i),p,'solution');
%     ff0 = obj.computefluxintegral(ZZ0,'solution');
%     R0 = (1/obj.tStep)*(obj.exactSolutionAll(:,i+1)-obj.exactSolutionAll(:,i))-ff1;
%     R1 = (1/obj.tStep)*(obj.exactSolutionAll(:,i+1)-obj.exactSolutionAll(:,i))-ff0;
%     else
%     ZZ1 = obj.unstructuredrecon(obj.exactSolutionAll(:,i+1),p,'error');
%     ff1 = obj.computefluxintegral(ZZ1,'solution');
%     ZZ0 = obj.unstructuredrecon(obj.exactSolutionAll(:,i),p,'error');
%     ff0 = obj.computefluxintegral(ZZ0,'solution');
%     
%     ZZ2 = obj.unstructuredrecon(obj.Uall(:,i),obj.pOrder,'solution');
%     ff2 = obj.computefluxintegral(ZZ2,'solution');
%     ZZ3 = obj.unstructuredrecon(obj.Uall(:,i+1),obj.pOrder,'solution');
%     ff3 = obj.computefluxintegral(ZZ3,'solution');
%     ZZ4 = obj.unstructuredrecon(obj.Uall(:,i),p,'error');
%     ff4 = obj.computefluxintegral(ZZ4,'solution');
%     ZZ5 = obj.unstructuredrecon(obj.Uall(:,i+1),p,'error');
%     ff5 = obj.computefluxintegral(ZZ5,'solution');
% %     r0 = R0;
% %     r1 = R1;
%     
%     
% %     R0 = (1/obj.tStep)*(obj.exactSolutionAll(:,i+1)-obj.exactSolutionAll(:,i))-ff1-ff3+ff5;
% %     R1 = (1/obj.tStep)*(obj.exactSolutionAll(:,i+1)-obj.exactSolutionAll(:,i))-ff0-ff2+ff4;
%   
%     DD0 = obj.residual(:,i)-ff5;
%     DD1 = obj.residual(:,i+1)-ff4;
% %     R1 = R1+zeros(size(R1))*h(2)^4;
% %     [r0+r1 R0+R1]
% % %     error('1')
% %     
%     end
    
    
    if(0 && time > obj.endTime/2)
       max(abs(r0+r1-R0-R1))
       max(abs(DD0+DD1+2*(obj.Uall(:,i+1)-obj.Uall(:,i))/k))
%            [r0+r1 R0+R1]
%        figure
%        plot(x,R0+R1);
obj.Uall(floor(obj.nCells/2),:)
       error('1')
       
    end
%     
%     Z = obj.unstructuredrecon(u,p,eqn);
%     f = obj.computefluxintegral(Z,eqn);
%     ZZt = obj.unstructuredrecon(u+obj.Uall(:,i),p,'solution');
%     ff = obj.computefluxintegral(ZZt,'solution');
%     ZZt = obj.unstructuredrecon(obj.Uall(:,i),p,'solution');
%     fff = obj.computefluxintegral(ZZt,'solution');
% 
%     v1 = obj.exactSolutionAll(:,1);
%     v2 = obj.exactSolutionAll(:,2);
%     Q = obj.unstructuredrecon(v1,p,'solution');
%     v3 = obj.computefluxintegral(Q,'solution');
%     Q = obj.unstructuredrecon(v2,p,'solution');
%     v4 = obj.computefluxintegral(Q,'solution');
%     Q = obj.unstructuredrecon(obj.Uall(:,2),p,'solution');
%     v5 = obj.computefluxintegral(Q,'solution');
    
%     [obj.Uall(:,1) obj.exactSolutionAll(:,1)]
%     [obj.Uall(:,1)+(k/2)*(v3+v5) obj.Uall(:,2)]
%     [v2-obj.Uall(:,2) (k/2)*(v4-v5)+(k/2)*((v2-v1)/k-v4+(v2-v1)/k-v5)]
%     error('1')
% if( max(abs(f-(ff-fff)))>1e-12 )
%     i
%         [f ff fff ff-fff u obj.Uall(:,i)]
%    error('1') 
% end

end
u0 = u;
d = 0;
g = 1;


    
% u* = u0-(1-(k/2)J)^-1*(u0-(un+(k/2)*(fn+f0)))
while(max(abs(g)) > 1e-12)
    
    obj.curTime = obj.curTime + obj.tStep;
    Z = obj.unstructuredrecon(u,p,eqn);
    f1 = obj.computefluxintegral(Z,eqn);
    obj.curTime = obj.curTime - obj.tStep;
    Z = obj.unstructuredrecon(u0,p,eqn);
    f0 = obj.computefluxintegral(Z,eqn);
    
%       if (strcmp(eqn,'error')==1 && max(abs(u)) > 1e-12 )
%           ZZ0 = obj.unstructuredrecon(u+obj.Uall(:,i+1),p,'solution');
%           ff0 = obj.computefluxintegral(ZZ0,'solution');
%           ZZ1 = obj.unstructuredrecon(obj.Uall(:,i+1),p,'solution');
%           ff1 = obj.computefluxintegral(ZZ1,'solution');
%         [u f0 f1 ff0-ff1 R0 R1]
% %         error('1')
%       end
    
    
    g = (u-u0)-0.5*k*(f0+f1);
    
    
    
    R = obj.computefluxjacobian(u,eqn);
    J = eye(N+2,N+2)-(k/2)*R;
    
    if (strcmp(eqn,'error')==1)
        
%         [u f0 f1]
%         error('1')

       g = g-0.5*k*(R0+R1); 
       
%         istep = round(obj.curTime/obj.tStep) +1;
%         K = obj.primalJacobian;        
%         KK = K(2:N+1,2:N+1);
%         g(2:N+1) =  g(2:N+1) - (1/12)*obj.tStep^3*KK^3*obj.Uall(2:N+1,istep);
        
%        R = obj.computefluxjacobian(u+obj.Uall(:,i),eqn);
%        RR = obj.computefluxjacobian(u+obj.Uall(:,i+1),eqn);
%        J = eye(N+2,N+2)-(k/2)*(R+RR);
% max(abs(g))
    end
    u(2:N+1) = u(2:N+1)-(J(2:N+1,2:N+1))\g(2:N+1);
%     [obj.curTime max(abs(g))]
    
end
% if(strcmp(eqn,'error')==1 && obj.curTime > 0.089)
% % [obj.exactSolutionAll(:,2)-obj.Uall(:,2) u obj.exactSolutionAll(:,2)-obj.Uall(:,2)-u]
% [obj.exactSolutionAll(:,11)-obj.Uall(:,11) u obj.exactSolutionAll(:,11)-obj.Uall(:,11)-u]
% error('1')
% end

uu = u;
uu(N+2) = NaN;

d = max(d,abs((uu-u0))/k);

if(strcmp(eqn,'error')==1) % check epsilon_1 for 244? e1 = e0 + (k/2)( f(u1+e1)-f(u1) ) + (k/2) (R0 + R1)
%      Z = obj.unstructuredrecon(u0,p,eqn);
%     g0 = obj.computefluxintegral(Z,eqn);
%     Z = obj.unstructuredrecon(uu,p,eqn);
%     obj.curTime = obj.curTime + obj.tStep;
%     g1 = obj.computefluxintegral(Z,eqn);
%     obj.curTime = obj.curTime - obj.tStep;
%    [u0 uu (k/2)*(g0+g1)+(k/2)*(R0+R1)]
   %g0 g1 R0 R1
%    if(i==1)
%       obj.higherprimalPI =  (k/2)*(g0+g1)+(k/2)*(R0+R1);
%    end
%    if(i==2)
%        e1 = obj.higherprimalPI;
%        [uu e1+(k/2)*(g0+g1)+(k/2)*(R0+R1)]
%           error('1')
%    end

end
    

if(strcmp(eqn,'solution')==1)
    
    i = round(obj.curTime/obj.tStep) +1;
    
    %     obj.Rall(:,2*i) = y;
    Z = obj.unstructuredrecon(u,p,eqn);
    f = obj.computefluxintegral(Z,eqn);
    obj.Rall(:,i+1) = f;
    
    J = obj.computefluxjacobian(u,eqn);
    %     obj.Uall(2:N+1,2*i) = J(2:N+1,2:N+1)\y(2:N+1);
%     obj.Uall(2:N+1,i+1) = J(2:N+1,2:N+1)\f(2:N+1);
    obj.Uall(2:N+1,i+1) = u(2:N+1);
    
    %     assert(obj.linearPhysics);
end
end

