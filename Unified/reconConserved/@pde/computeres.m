function [ R ] = computeres(obj,u,time,r)%u,x,h,N,f,r,phys,time,gsp )
%COMPUTERES Summary of this function goes here
%   Detailed explanation goes here
global TEND
if(abs(time-round(time))<1e-10)
   time = round(time); 
elseif(abs(time-TEND) < 1e-10)
   time = TEND;
end


% r = obj.rOrder;
N = obj.nCells;
h = obj.cellWidths;
phys = obj.physics;
f = obj.source;

eqn = 'residual';
if(r==obj.pOrder)
   eqn = 'solution'; 
end
Z = obj.unstructuredrecon(u,r,eqn);
%  if (time > 0)
%     Z
%     eqn
% 
%     error('1')
%  end



% %  [Z]=unstructuredrecon(u,x,h,N,NaN,NaN,r);
 R = zeros(N+2,1);
 %uxx = zeros(N+2,1);
 
% % % %   [left,right] = computeflux(Z,h,N,r,phys,'residual',obj);
  
 R= obj.computefluxintegral(Z,'residual');  
 
 



%   error('1')
 for i = 2:N+1

% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
  
% % %  [left,right] = computeflux(Z,h,i,N,r,phys);
%  sp = gsp(i);
% ut = fnval(fnder(sp),time);

global dUdt
global KK


ut = dUdt(i,round(time/KK)+1);

% if (time ==0)
% round(time/KK)+1
% end

% % % % R(i)= -ut+(right(i)-left(i))/h(i)-f(i);
R(i) = R(i) - ut;


if(strcmp(phys,'Burgers')==1)
   R(i) = -ut-(right(i)^2-left(i)^2)/(2*h(i))-f(i);  
end


 end
 
 
%   dUdt
%   error('1')
 
% switch r
% 
%         case 2    
%     [err,Z]=unstructuredrecon1(u,x,h,N,NaN,NaN);
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% for i = 2:N+1
%    
% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
%  
% end
% 
%     case 3    
%     [err,Z]=unstructuredrecon2(u,x,h,N,NaN,NaN);
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% for i = 2:N+1
%    
% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
%  
% end
% case 4
%     
%     
%     [err,Z]=unstructuredrecon3(u,x,h,N,NaN,NaN);
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% for i = 2:N+1
%   
%   
% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
% 
% end
% 
% 
% 
% case 5
%     [err,Z]=unstructuredrecon4(u,x,h,N,NaN,NaN);
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% 
% for i = 2:N+1
% 
% 
% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
% 
% 
% end
%     
% 
% case 6
%     [err,Z]=unstructuredrecon5(u,x,h,N,NaN,NaN);
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% 
%  for i = 2:N+1
% 
% [ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r,phys,uder,j,time,gsp);%%%
% 
% end
%     
%     otherwise 
%         r
%                assert(0==1)
% end
% 
% end

