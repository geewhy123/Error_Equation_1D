function [ J ] = computefluxjacobian( obj,u,x,h,N,p)
%COMPUTEFLUXJACOBIAN Summary of this function goes here
%   Detailed explanation goes here

J = zeros(N+2,N+2);
Z = obj.unstructuredrecon(u);
R = obj.computefluxintegral(Z);%,x,h,N,p)
I = eye(N);
ep = 1e-10;
u1 = NaN*ones(N+2,1);
for i = 2:N+1
%    for j = 2:N+1 
       u1(2:N+1,1) = u(2:N+1) + ep*I(:,i-1);
      
       Z1 = obj.unstructuredrecon(u1);
       R1=obj.computefluxintegral(Z1);%u1,x,h,N,p);
       R0=obj.computefluxintegral(Z);%u,x,h,N,p);
   J(2:N+1,i) = (R1(2:N+1)-R0(2:N+1))/ep;  
    
%    end
end



end

