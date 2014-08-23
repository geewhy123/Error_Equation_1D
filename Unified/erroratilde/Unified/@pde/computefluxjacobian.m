function [ J ] = computefluxjacobian( obj,u,eqn)
%COMPUTEFLUXJACOBIAN Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(eqn,'solution')==1)
    order = obj.pOrder;
elseif(strcmp(eqn,'error')==1)
    order = obj.qOrder;
else
assert(0);    
end


N = obj.nCells;
J = zeros(N+2,N+2);

Z = obj.unstructuredrecon(u,order,eqn);

       R0=obj.computefluxintegral(Z,eqn);%u,x,h,N,p);
% R = obj.computefluxintegral(Z);%,x,h,N,p)




I = eye(N);
ep = 1e1;
if(strcmp(obj.physics,'Poisson') ~= 1 &&strcmp(obj.physics,'Advection') ~= 1 )
%    fprintf('not linear problem, pick smaller epsilon for FD Jacobian')
   ep = 1e-8;
%    assert(0);
end

u1 = NaN*ones(N+2,1);

for i = 2:N+1
%    for j = 2:N+1 
       u1(2:N+1,1) = u(2:N+1) + ep*I(:,i-1);

 
       
       Z1 = obj.unstructuredrecon(u1,order,eqn);
 
       R1=obj.computefluxintegral(Z1,eqn);%u1,x,h,N,p);
 
%        R0=obj.computefluxintegral(Z,eqn);%u,x,h,N,p);

       
%     if(strcmp(eqn,'error')==1)
%     R0
% error('1')
%  end

       
       
%        R1-R0
%        error('1')
   J(2:N+1,i) = (R1(2:N+1)-R0(2:N+1))/ep;  
    
%    end
end




end

