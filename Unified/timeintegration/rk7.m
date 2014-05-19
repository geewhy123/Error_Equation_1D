function [uu,d] = rk7(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR)
%RK1 Summary of this function goes here
%   Detailed explanation goes here

% b = [7/90 0 32/90 12/90 32/90 7/90];
% A=[

%time = time-k;
c = [0 1/6 1/3 1/2 2/11 2/3 6/7 0 1 ];
Zu=  NaN*ones(p,N+2,9);
Ubar = NaN*ones(N+2,9);

val = NaN*ones(N+2,9);
if (strcmp(eqn,'error')==1)


timesbet = c*k+time;
   
global TEND
% if(abs(timesbet(end)-round(time))<1e-10)
%    timesbet(end) = round(timesbet(end)); 
% else
    if(abs(timesbet(end)-TEND) < 1e-10)
   timesbet(end) = TEND;
end

for kk = 2:N+1
val(kk,1:length(timesbet)) = fnval(timesbet,Rsp(kk));
end
%val

end

%nonlinear error
  if (strcmp(phys,'Burgers')==1)
% if((nargin < 12) || (isnan(time))|| isnan(j))%primal and error step
     if(~isnan(time))
         Ubar = zeros(N+2,9);
         global UU
         global M        
             T=(0:1:M)*k;
             
         for i = 2:N+1
            Usp(i) = spapi(6,T,UU(i,:));
            
%             for steps = 1:9
            Ubar(i,1:9) = fnval(Usp(i),time+c*k);
%             end
         
         end
         
         for steps = 1:9
            [Zu(:,:,steps)]=unstructuredrecon(Ubar(:,steps),x,h,N,NaN,NaN,p); 
         end
        
     end
%   end
  end


% phi = zeros(N+2,1);
% phiII = zeros(N+2,1);
% phiIII = zeros(N+2,1);
% phiIV = zeros(N+2,1);
% phiV = zeros(N+2,1);
% phiVI = zeros(N+2,1);
% phiVII = zeros(N+2,1);
% phiVIII = zeros(N+2,1);
% phiIX = zeros(N+2,1);
% uII = zeros(N+2,1);
% uIII = zeros(N+2,1);
% uIV = zeros(N+2,1);
% uV = zeros(N+2,1);
% uVI = zeros(N+2,1);
% uVII = zeros(N+2,1);
% uVIII = zeros(N+2,1);
% uIX = zeros(N+2,1);
unew = u;%zeros(N+2,1);
 uu = zeros(N+2,1);

 A = [1/6 0 0 0 0 0 0 0 0;
     0 1/3 0 0 0 0 0 0 0;
     1/8 0 3/8 0 0 0 0 0 0;
     148/1331 0 150/1331 -56/1331 0 0 0 0 0;
     -404/243 0 -170/27 4024/1701 10648/1701 0 0 0 0;
     2466/2401 0 1242/343 -19176/16807 -51909/16807 1053/2401 0 0 0;
     5/154 0 0 96/539 -1815/20384 -405/2464 49/1144 0 0;
     -113/32 0 -195/22 32/7 29403/3584 -729/512 1029/1408 21/16 0;
     0 0 0 32/105 1771561/6289920 243/2560 16807/74880 77/1440 11/270] ;
 
 phi = zeros(N+2,length(c));
 for steps = 1:length(c);
    Z = unstructuredrecon(unew,x,h,N,BCLeft,uL,BCRight,uR,p);
    phi(:,steps) = timestep(eqn,Z,f,k,h,N,p,phys,time+c(steps)*k,Rsp,Zu(:,:,steps),val(:,steps));
    
    unew = u+phi*A(steps,:)'*k;
    
 end

 d = 0;
 uu = unew;
 d = max(d,abs((uu-u))/k); 
%  
%  
% [Z]=unstructuredrecon(u,x,h,N,NaN,NaN,p);
% d = 0;
% % for i = 2:N+1
% phi = timestep(eqn,Z,f,k,h,N,p,phys,time+c(1)*k,Rsp,Zu(:,:,1),val(:,1));
% unew = u+(k/6)*phi;
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiII = timestep(eqn,Z,f,k,h,N,p,phys,time+c(2)*k,Rsp,Zu(:,:,2),val(:,2));
% unew = u+(k/3)*(phiII);
%    
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiIII = timestep(eqn,Z,f,k,h,N,p,phys,time+c(3)*k,Rsp,Zu(:,:,3),val(:,3));
% unew = u+(k/8)*(phi+3*phiIII);
%    
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiIV = timestep(eqn,Z,f,k,h,N,p,phys,time+c(4)*k,Rsp,Zu(:,:,4),val(:,4));
% unew = u+(k/1331)*(148*phi+150*phiIII-56*phiIV);
%    
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiV = timestep(eqn,Z,f,k,h,N,p,phys,time+c(5)*k,Rsp,Zu(:,:,5),val(:,5));
% unew = u+(k/1701)*(-2828*phi-10710*phiIII+4024*phiIV+10648*phiV);
%    
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiVI = timestep(eqn,Z,f,k,h,N,p,phys,time+c(6)*k,Rsp,Zu(:,:,6),val(:,6));
% unew = u+(k/16807)*(17262*phi+60858*phiIII-19176*phiIV-51909*phiV+7371*phiVI);
%    
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiVII = timestep(eqn,Z,f,k,h,N,p,phys,time+c(7)*k,Rsp,Zu(:,:,7),val(:,7));
% unew = u+(k)*((5/154)*phi+(96/539)*phiIV-(1815/20384)*phiV-(405/2464)*phiVI+(49/1144)*phiVII);
%    
% %d = max(d,abs(delt)); 
% % end
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiVIII = timestep(eqn,Z,f,k,h,N,p,phys,time+c(8)*k,Rsp,Zu(:,:,8),val(:,8));
% unew = u+(k)*((-113/32)*phi-(195/22)*phiIII+(32/7)*phiIV+(29403/3584)*phiV-(729/512)*phiVI+(1029/1408)*phiVII+(21/16)*phiVIII);
%    
% %d = max(d,abs(delt)); 
% %  end
% 
% 
% [Z]=unstructuredrecon(unew,x,h,N,NaN,NaN,p);
% % for i = 2:N+1        
%        
% phiIX = timestep(eqn,Z,f,k,h,N,p,phys,time+c(9)*k,Rsp,Zu(:,:,9),val(:,9));
% 
% uu = u+(k)*((32/105)*phiIV+(1771561/6289920)*phiV+(243/2560)*phiVI+(16807/74880)*phiVII+(77/1440)*phiVIII+(11/270)*phiIX);   
%d = max(d,abs(((32/105)*phiIV+(1771561/6289920)*phiV+(243/2560)*phiVI+(16807/74880)*phiVII+(77/1440)*phiVIII+(11/270)*phiIX))); 
% end
uu(N+2) = NaN;
   


end

