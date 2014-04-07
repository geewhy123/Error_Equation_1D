function [ upr,upl,phi] = reconfluxsoln( u,Z,f,k,h,N,p,phys,uder,j,time,gsp,Rsp,Zu)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here
phi = zeros(N+2,1);
phi(1) = NaN;
phi(N+2) = NaN;


%  if((nargin < 12) || (isnan(time))|| isnan(j))
%         if(~isnan(time))% error eqn
% Rsp(3)
% val(2:N+1) = fnval(time,Rsp(2:N+1));
%         end
%  end

for i = 2:N+1


global xx
global PHI
global TEND
if(abs(time-round(time))<1e-10)
   time = round(time); 
elseif(abs(time-TEND) < 1e-10)
   time = TEND;
end

global UU
global M

[left,right] = computeflux(Z,h,i,N,p,phys);
ur1 = right;
ur2 = right;
upr = right;
ul1 = left;
ul2 = left;
upl = left;



if(strcmp(phys,'Poisson')==1)
   
        phi(i)= (upr-upl)/h(i)-f(i);%Poisson
  
elseif(strcmp(phys,'Advection')==1)
    

        
        phi(i)= (ur2-ul1)/h(i)-f(i); % primal
   

    
elseif(strcmp(phys,'Burgers')==1)
    

        
             phi(i) = -(ur1^2-ul2^2)/(2*h(i))-f(i);%burgers
    
end


end

end
