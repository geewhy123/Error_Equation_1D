function [uu,d] = rk7(eqn,u,x,f,k,h,N,p,t,phys,uder,j,time,gsp,Rsp)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

% b = [7/90 0 32/90 12/90 32/90 7/90];
% A=[

%time = time-k;
c = [0 1/6 1/3 1/2 2/11 2/3 6/7 0 1 ];
Zu=  NaN*ones(p,N+2,9);
Ubar = NaN*ones(N+2,9);
%nonlinear error
  if (strcmp(phys,'Burgers')==1)
  if((nargin < 12) || (isnan(time))|| isnan(j))%primal and error step
     if(~isnan(time))
         Ubar = zeros(N+2,9);
         global UU
         global M        
             T=(0:1:M)*k;
             
         for i = 2:N+1
            Usp(i) = spapi(6,T,UU(i,:));
            
            for steps = 1:9
            Ubar(i,steps) = fnval(Usp(i),time+c(steps)*k);
            end
         
         end
         
         for steps = 1:9
            [Zu(:,:,steps)]=unstructuredrecon(Ubar(:,steps),x,h,N,NaN,NaN,p); 
         end
        
     end
  end
  end
  
 

phi = zeros(N+2,1);
phiII = zeros(N+2,1);
phiIII = zeros(N+2,1);
phiIV = zeros(N+2,1);
phiV = zeros(N+2,1);
phiVI = zeros(N+2,1);
phiVII = zeros(N+2,1);
phiVIII = zeros(N+2,1);
phiIX = zeros(N+2,1);
uII = zeros(N+2,1);
uIII = zeros(N+2,1);
uIV = zeros(N+2,1);
uV = zeros(N+2,1);
uVI = zeros(N+2,1);
uVII = zeros(N+2,1);
uVIII = zeros(N+2,1);
uIX = zeros(N+2,1);
uu = zeros(N+2,1);

[Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
d = 0;
% for i = 2:N+1
       
[upr,upl,phi] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(1)*k,gsp,Rsp,Zu(:,:,1));
uII = u+(k/6)*phi;
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiII] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(2)*k,gsp,Rsp,Zu(:,:,2));
uIII = u+(k/3)*(phiII);
   
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uIII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiIII] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(3)*k,gsp,Rsp,Zu(:,:,3));
uIV = u+(k/8)*(phi+3*phiIII);
   
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uIV,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiIV] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(4)*k,gsp,Rsp,Zu(:,:,4));
uV = u+(k/1331)*(148*phi+150*phiIII-56*phiIV);
   
% end

[Z]=unstructuredrecon(uV,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiV] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(5)*k,gsp,Rsp,Zu(:,:,5));
uVI = u+(k/1701)*(-2828*phi-10710*phiIII+4024*phiIV+10648*phiV);
   
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uVI,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiVI] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(6)*k,gsp,Rsp,Zu(:,:,6));
uVII = u+(k/16807)*(17262*phi+60858*phiIII-19176*phiIV-51909*phiV+7371*phiVI);
   
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uVII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiVII] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(7)*k,gsp,Rsp,Zu(:,:,7));
uVIII = u+(k)*((5/154)*phi+(96/539)*phiIV-(1815/20384)*phiV-(405/2464)*phiVI+(49/1144)*phiVII);
   
%d = max(d,abs(delt)); 
% end

[Z]=unstructuredrecon(uVIII,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiVIII] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(8)*k,gsp,Rsp,Zu(:,:,8));
uIX = u+(k)*((-113/32)*phi-(195/22)*phiIII+(32/7)*phiIV+(29403/3584)*phiV-(729/512)*phiVI+(1029/1408)*phiVII+(21/16)*phiVIII);
   
%d = max(d,abs(delt)); 
%  end


[Z]=unstructuredrecon(uIX,x,h,N,NaN,NaN,p);
% for i = 2:N+1        
       
[upr,upl,phiIX] = timestep(eqn,u,Z,f,k,h,N,p,phys,uder,j,time+c(9)*k,gsp,Rsp,Zu(:,:,9));

uu = u+(k)*((32/105)*phiIV+(1771561/6289920)*phiV+(243/2560)*phiVI+(16807/74880)*phiVII+(77/1440)*phiVIII+(11/270)*phiIX);   
d = max(d,abs(((32/105)*phiIV+(1771561/6289920)*phiV+(243/2560)*phiVI+(16807/74880)*phiVII+(77/1440)*phiVIII+(11/270)*phiIX))); 
% end
uu(N+2) = NaN;
   


end

