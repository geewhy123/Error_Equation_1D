function [phi] = reconfluxerror( Z,f,k,h,N,p,phys,time,Rsp,Zu,val)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here
phi = zeros(N+2,1);
phi(1) = NaN;
phi(N+2) = NaN;




for i = 2:N+1


global TEND
if(abs(time-round(time))<1e-10)
   time = round(time); 
elseif(abs(time-TEND) < 1e-10)
   time = TEND;
end



[left,right] = computeflux(Z,h,i,N,p,phys);
ur1 = right;
ur2 = right;
upr = right;
ul1 = left;
ul2 = left;
upl = left;



if(strcmp(phys,'Poisson')==1)
   
            %time
            
%            sp = Rsp(i);
     
            f(i) = -1*val(i);%getRes(time,k,i,sp);
    
            %val(i)
            %getRes(time,k,i,sp)
           
            %if(abs(val(i) - getRes(time,k,i,sp))>1e-3)
             %   val
             %   fnval(sp,time)
             %   error('1');
            %end
            
            
        phi(i)= (upr-upl)/h(i)-f(i);%Poisson

    
elseif(strcmp(phys,'Advection')==1)
    
    
            
%             sp = Rsp(i);
            f(i) = -1*val(i);%getRes(time,k,i,sp);
        
        
        phi(i)= (ur2-ul1)/h(i)-f(i); % primal
        

    
elseif(strcmp(phys,'Burgers')==1)
    %nonlinear error

% if((nargin < 12) || (isnan(time))|| isnan(j))%primal and error step
%    %  if(~isnan(time))
%          Ubar = zeros(N+2,1);
%          global UU
%          global M        
%              T=(0:1:M)*k;
%              
% %          for i = 2:N+1
%             Usp(i) = spapi(6,T,UU(i,:));
%             
% %             for steps = 1:9
% %             Ubar(i,steps) = fnval(Usp(i),time+c(steps)*k);
% %             end
%          Ubar(i) = fnval(Usp(i),time);
% %          end
%          
% %          for steps = 1:9
%             [Zu(:,:,steps)]=unstructuredrecon(Ubar(:,steps),x,h,N,NaN,NaN,p); 
%          end
%         
%     % end
% %   end
% %   end
  
  
     
            %time
%             sp = Rsp(i);
         
            
            f(i) = -1*val(i);%getRes(time,k,i,sp);
            

%             switch i
%                case 2
%                 Y = Zu(:,i);
%                 Yr = Zu(:,i+1);
%                 Yl = Zu(:,N+1);
%                case N+1
%                 Y = Zu(:,i);
%                 Yr = Zu(:,2);
%                 Yl = Zu(:,i-1);
%                otherwise
%                 Y = Zu(:,i);
%                 Yr = Zu(:,i+1);
%                 Yl = Zu(:,i-1);
%                end
%             if (p==2)
%                
%                   uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2);
%                   uhr1 = Y(1) + Y(2)*(h(i)/2);
%         
%                   uhl1 = Y(1)+Y(2)*(-h(i)/2);
%                   uhl2 = Yl(1) + Yl(2)*(h(i-1)/2);
%                   
%             elseif(p==4)
%                 uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3;
%                 uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3;
%                 uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3;
%                 uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3;
%             elseif(p==5)
%                 uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3+Yr(5)*(-h(i+1)/2)^4;
%                 uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3 + Y(5)*(h(i)/2)^4;
%                 uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3+Y(5)*(-h(i)/2)^4;
%                 uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3 +Yl(5)*(h(i-1)/2)^4; 
%             elseif(p==6)
%                 uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3+Yr(5)*(-h(i+1)/2)^4 +Yr(6)*(-h(i+1)/2)^5;
%                 uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3 + Y(5)*(h(i)/2)^4 +Y(6)*(h(i)/2)^5;
%                 uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3+Y(5)*(-h(i)/2)^4 +Y(6)*(-h(i)/2)^5;
%                 uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3 +Yl(5)*(h(i-1)/2)^4 + Yl(6)*(h(i-1)/2)^5; 
%             else
%                 error('2')
%             end

[left,right] = computeflux(Zu,h,i,N,p,phys);

%             f(i) = f(i) + (ur1*(uhr1)-ul2*(uhl2))/h(i);
f(i) = f(i) + (ur1*(right)-ul2*(left))/h(i);
      
        
             phi(i) = -(ur1^2-ul2^2)/(2*h(i))-f(i);%burgers
                    

 
   
end


end

end
