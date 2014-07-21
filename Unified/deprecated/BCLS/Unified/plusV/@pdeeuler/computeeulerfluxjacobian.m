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
%%%%%%
%  m = obj.hOrder;
% % [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
%   [Z3] = higherunstructuredreconeuler (obj,v(:,3),m,'solution',3);                          
%   [Z1] = higherunstructuredreconeuler (obj,v(:,1),m,'solution',1); 
%   [Z2] = higherunstructuredreconeuler (obj,v(:,2),m,'solution',2);
%                Z = [Z1; Z2;Z3];
% 
% 
% %%%%%%


  [phi1,phi2,phi3]=obj.computeeulerfluxintegral(Z,eqn);


R0 = NaN*ones(3*N+2,1);
R1 = NaN*ones(3*N+2,1);
 R0(2:3:3*N-1) = phi1(2:N+1);
 R0(3:3:3*N) = phi2(2:N+1);
 R0(4:3:3*N+1) = phi3(2:N+1);

u = NaN*ones(size(v));

if(strcmp(eqn,'error')==1 && 0)
Upe =zeros(N+2,3);
Vpe = zeros(N+2,3);
% %     for j = 2:N+1
% %     u(j,1) = v(j,1);
% %     u(j,2) = v(j,2);
% %     u(j,3) = v(j,3);
% %     end
% % incomplete here
    Vpe = v+obj.convSolutionV;
    for j = 2:N+1
    [Upe(j,1),Upe(j,2),Upe(j,3)] = toconservedvars(Vpe(j,1),Vpe(j,2),Vpe(j,3));
    end
    u = Upe-obj.convSoln
v
% error('1')

%     Vpe = u+V;
else
for j = 2:N+1
    [u(j,1),u(j,2),u(j,3)] = toconservedvars(v(j,1),v(j,2),v(j,3));
end
end

U = NaN*ones(3*N+2,1);
U(2:3:3*N-1) = u(2:N+1,1);
U(3:3:3*N) = u(2:N+1,2);
U(4:3:3*N+1) = u(2:N+1,3);


% U
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

if(strcmp(eqn,'error')==1 && 0)

% %     for j = 2:N+1
% %     v1(j,1) = u1(j,1);
% %     v1(j,2) = u1(j,2);
% %     v1(j,3) = u1(j,3);
% %     end
    Upe = u1+obj.convSoln;
    for j = 2:N+1
    [Vpe(j,1),Vpe(j,2),Vpe(j,3)] = toprimitivevars(Upe(j,1),Upe(j,2),Upe(j,3));
    end
    v1 = Vpe-obj.convSolutionV;
    

% incomplete here
else
for j = 2:N+1
[v1(j,1),v1(j,2),v1(j,3)] = toprimitivevars(u1(j,1),u1(j,2),u1(j,3));
end
end

       Z1 = obj.unstructuredrecon(v1,order,eqn);

% %%%%%%
% % [Z] = obj.unstructuredrecon(V,p,'solution');%u,x,h,N,NaN,NaN,p);
%   [Z3] = higherunstructuredreconeuler (obj,v1(:,3),m,'solution',3);                          
%   [Z0] = higherunstructuredreconeuler (obj,v1(:,1),m,'solution',1); 
%   [Z2] = higherunstructuredreconeuler (obj,v1(:,2),m,'solution',2);
%                Z1 = [Z0; Z2;Z3];
% 
% 
% %%%%%%


 
       [phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z1,eqn);%u1,x,h,N,p);


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

           [R1 R0 R1-R0];
% if(strcmp(eqn,'error')==1)
%         error('1')
% end
%        end
   J(2:3*N+1,i) = (R1(2:3*N+1)-R0(2:3*N+1))/ep;  

   
for k = 2:3*N+1
  if(abs(R1(k)-R0(k)) < 1e-13)
    J(k,i) = 0;   
  end
end

%     if( i == 2 && strcmp(eqn,'error')==1)
%      [U U1]
%         [R1(2:3*N+1) R0(2:3*N+1)]
%      error('1')
%     end
%    end
end



end

