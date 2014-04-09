function [phi] = reconfluxsoln( Z,f,h,N,p,phys,time)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here
phi = zeros(N+2,1);
phi(1) = NaN;
phi(N+2) = NaN;


[left,right] = computeflux(Z,h,N,p,phys);

% for i = 2:N+1


% global TEND
% if(abs(time-round(time))<1e-10)
%    time = round(time); 
% elseif(abs(time-TEND) < 1e-10)
%    time = TEND;
% end


% 
% [left,right] = computeflux(Z,h,i,N,p,phys);
ur1 = right;
ur2 = right;
upr = right;
ul1 = left;
ul2 = left;
upl = left;



if(strcmp(phys,'Poisson')==1)
   
        phi= (upr-upl)./h-f;%Poisson
  
elseif(strcmp(phys,'Advection')==1)
    

        
        phi= (ur2-ul1)./h-f; % primal
   

    
elseif(strcmp(phys,'Burgers')==1)
    

        
             phi = -(ur1.^2-ul2.^2)./(2*h)-f;%burgers
    
end


end

% end
