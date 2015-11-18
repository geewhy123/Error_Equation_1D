function [uu,d] = rk7(eqn,u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj)
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
%    error('1')
end

for kk = 2:N+1
val(kk,1:length(timesbet)) = fnval(timesbet,obj.Rsp(kk));
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
% %             [Zu(:,:,steps)]=unstructuredrecon(Ubar(:,steps),x,h,N,NaN,NaN,p); 
            Zu(:,:,steps) = obj.unstructuredrecon(Ubar(:,steps),p);
         end
        
     end

%   end
  end

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
 for steps = 1:length(c)
%     Z = unstructuredrecon(unew,x,h,N,BCLeft,uL,BCRight,uR,p);
  
Z = obj.unstructuredrecon(unew,p,eqn);

% Z
% obj.reconplot(Z)
% error('1')
     
% % % %     phi(:,steps) = timestep(eqn,Z,f,k,h,N,p,phys,time+c(steps)*k,Rsp,Zu(:,:,steps),val(:,steps));
% if(steps==length(c))
% %     A
% %     k
% %     phi
% %     unew
%    Z
%    error('1')
% end
if(strcmp(eqn,'solution')==1)
    obj.curTime = obj.curTime + (steps-1)*k/length(c);
phi(:,steps) = obj.computefluxintegral(Z,eqn);

%  phi(:,steps)
%  error('1')



elseif(strcmp(eqn,'error')==1)
        obj.errorSource = -1*val(:,steps);
        
%         obj.curTime = obj.curTime + (steps-1)*k/length(c);
         phi(:,steps) = obj.computefluxintegral(Z,'error');
%         phi
%         obj.errorSource
%         error('1')
%       phi(:,steps) = timestep(eqn,Z,f,k,h,N,p,phys,time+c(steps)*k,Rsp,Zu(:,:,steps),val(:,steps),obj);

else
    error('2')
end
%     Z
%     phi
%     error('1')

% error('1')
unew = u+phi*A(steps,:)'*k;
%     (phi)
%     (A(steps,:)')
%     phi*A(steps,:)'*k
%     error('1')
U(steps,:) = phi*A(steps,:)';    
 end
 
% if(strcmp(eqn,'error')==1)
% % [u unew]
% % U'
% % phi
% % A(steps,:)'
% % k
%    error('1') 
% end
 d = 0;
 uu = unew;
 d = max(d,abs((uu-u))/k); 

uu(N+2) = NaN;
   


end

