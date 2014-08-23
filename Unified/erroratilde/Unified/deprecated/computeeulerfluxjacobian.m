function [ J ] = computeeulerfluxjacobian( obj,v,eqn)
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
J = zeros(3*N+2,3*N+2);

Z = obj.unstructuredrecon(v,order,eqn);
% Z
% error('1')

%        R0=obj.computefluxintegral(Z,eqn);%u,x,h,N,p);
% R = obj.computefluxintegral(Z);%,x,h,N,p)

[phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,eqn);

R0 = NaN*ones(3*N+2,1);
R1 = NaN*ones(3*N+2,1);
R0(2:3:3*N-1) = phi1(2:N+1);
R0(3:3:3*N) = phi2(2:N+1);
R0(4:3:3*N+1) = phi3(2:N+1);

u = NaN*ones(size(v));
for j = 2:N+1
[u(j,1),u(j,2),u(j,3)] = toconservedvars(v(j,1),v(j,2),v(j,3));
end

U = NaN*ones(3*N+2,1);
U(2:3:3*N-1) = u(2:N+1,1);
U(3:3:3*N) = u(2:N+1,2);
U(4:3:3*N+1) = u(2:N+1,3);
% error('1')



I = eye(3*N);
ep = 1e0;
if(strcmp(obj.physics,'Poisson') ~= 1 &&strcmp(obj.physics,'Advection') ~= 1 )
%    fprintf('not linear problem, pick smaller epsilon for FD Jacobian')
   ep = 1e-8;
%    assert(0);
end

U1 = NaN*ones(3*N+2,1);
u1 = NaN*ones(N+2,3);
v1 = u1;
for i = 2:3*N+1
%    for j = 2:N+1 

       U1(2:3*N+1,1) = U(2:3*N+1) + ep*I(:,i-1);

%  error('1')
       
u1(2:N+1,1) = U1(2:3:3*N-1) ;
u1(2:N+1,2) = U1(3:3:3*N) ;
u1(2:N+1,3) = U1(4:3:3*N+1);

for j = 2:N+1
[v1(j,1),v1(j,2),v1(j,3)] = toprimitivevars(u1(j,1),u1(j,2),u1(j,3));
end

       Z1 = obj.unstructuredrecon(v1,order,eqn);
 
       [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z1,eqn);%u1,x,h,N,p);

     R1(2:3:3*N-1) = phi1(2:N+1);
     R1(3:3:3*N) = phi2(2:N+1);
     R1(4:3:3*N+1) = phi3(2:N+1);
     
%        R0=obj.computefluxintegral(Z,eqn);%u,x,h,N,p);

       
%     if(strcmp(eqn,'error')==1)
%     R0
% error('1')
%  end

       
%        if(i==4)
%            Z1-Z
%            U1-U
%          [R1-R0]
%         error('1')
%        end
   J(2:3*N+1,i) = (R1(2:3*N+1)-R0(2:3*N+1))/ep;  
    
%    end
end




end

