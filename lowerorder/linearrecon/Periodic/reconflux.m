function [ upr,upl,phi] = reconflux( u,Z,f,k,h,i,N,p,phys,uder,j,time,gsp,Rsp,Zu)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here
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

switch i
    case 2
        y = Z(:,i);
        yr = Z(:,i+1);
        yl = Z(:,N+1);
    case N+1
        y = Z(:,i);
        yr = Z(:,2);
        yl = Z(:,i-1);
    otherwise
        y = Z(:,i);
        yr = Z(:,i+1);
        yl = Z(:,i-1);
    
end    

switch p
    case 2
        upr1 = y(2);
        upr2 = yr(2);
        ur2 = yr(1)+yr(2)*(-h(i+1)/2);
        ur1 = y(1) + y(2)*(h(i)/2);
        upr = (upr1+upr2)/2 + (.2/((h(i)+h(i+1))/2))*(ur2-ur1);
        
          
        upl1 = y(2);
        upl2 = yl(2);
        ul1 = y(1)+y(2)*(-h(i)/2);
        ul2 = yl(1) + yl(2)*(h(i-1)/2);
        upl = (upl1+upl2)/2 + (.2/((h(i)+h(i-1))/2))*(ul1-ul2);


    case 3
        upr1 = y(2)+2*y(3)*h(i)/2; % u_i+1/2 using recon in i 
        upr2 = yr(2)+2*yr(3)*(-h(i+1)/2);
        ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2;
        ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2;


        upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur2-ur1);
        upl1 = y(2)+2*y(3)*(-h(i)/2);
        upl2 = yl(2)+2*yl(3)*h(i-1)/2;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 ;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 ;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
        

    case 4
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
        upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
        ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3;
        ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3;
        upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur2-ur1);

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;

        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ul1-ul2);
             

    case 5
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3;
        upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
        ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4;
        ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3 + y(5)*(h(i)/2)^4;
        upr = (upr1+upr2)/2;

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4; 
        upl = (upl1+upl2)/2;


    case 6
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3+5*y(6)*(h(i)/2)^4;
        upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3 + 5*yr(6)*(-h(i+1)/2)^4;
        ur2 = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4 +yr(6)*(-h(i+1)/2)^5;
        ur1 = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3 + y(5)*(h(i)/2)^4 +y(6)*(h(i)/2)^5;
        upr = (upr1+upr2)/2;

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4; 
        ul1 = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4 +y(6)*(-h(i)/2)^5;
        ul2 = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3 +yl(5)*(h(i-1)/2)^4 + yl(6)*(h(i-1)/2)^5; 
        upl = (upl1+upl2)/2;
           
        
    otherwise
        error('1')
end


if(strcmp(phys,'Poisson')==1)
    if(nargin < 12 || (isnan(time))||isnan(j))
        if(~isnan(time))
     
            time
            
            sp = Rsp(i);
     
            f(i) = -1*getRes(time,k,i,sp);
        end
        phi= (upr-upl)/h(i)-f(i);%Poisson
    else
       %%% phi= -uder(i,j)+(upr-upl)/h(i)-f(i);%Poisson
       
sp = gsp(i);
ut = fnval(fnder(sp),time);
phi= -ut+(upr-upl)/h(i)-f(i);%Poisson
    end   
    
elseif(strcmp(phys,'Advection')==1)
    
    
    if((nargin < 12) || (isnan(time))|| isnan(j))
        if(~isnan(time))% error eqn
            
            sp = Rsp(i);
            f(i) = -1*getRes(time,k,i,sp);
        end
        
        phi= (ur2-ul1)/h(i)-f(i); % primal
        
    else%res

sp = gsp(i);
ut = fnval(fnder(sp),time);

phi= -ut+(ur2-ul1)/h(i)-f(i);

   if(abs(time-0.04)<1e-3)
       ut
       ur2
       ul1
       -ut+(ur2-ul1)/h(i)

           PHI(i) = phi;
           PHI(i)
           
           if(i==N+1)
               PHI(N+2) = NaN;
               PHI
               max(abs(PHI))
              %plot(xx(2:N+1),PHI(2:N+1))
              %clear global
              
              %error('1')
           end

        end

% if(norm(phi) >1)
%    time
%    ut
%    error('3')
% end


% if(~isnan(j) && abs(uder(i,j)-ut)/ut > 1e-4)
%    
%    uder(i,j)
%    ut
%    i
%    j
%    time-1
%    fnval(fnder(sp),1+1e-8)
%   error('1')
% end
%     
% % end
%         if((i==11) && abs(time-0.2)<1e-3)
%             time
%             j = round(time/k+1)
%             uder(i,j)
%             ut
%             uder(i,j)-ut
%             (ur2-ul1)/h(i)
%             f(i)
%            error('1');
%             phi;
%         end
    end    
    
    
elseif(strcmp(phys,'Burgers')==1)
    
    if((nargin < 12) || (isnan(time))|| isnan(j))%primal and error step
        if(~isnan(time))%error step
     
            %time
            sp = Rsp(i);
         
            
            f(i) = -1*getRes(time,k,i,sp);
            
%             T=(0:1:M)*k;
%             Usp1 = spapi(6,T,UU(i,:));
%             
%             if (i==N+1)
%             Usp2 = spapi(6,T,UU(2,:));
%             else
%             Usp2 = spapi(6,T,UU(i+1,:));
%             end
%             uhr = fnval(Usp2,time);
%             uhl = fnval(Usp1,time);
  %          f(i) = f(i) -(ur1*(uhr)-ur2*(uhl))/h(i);
            
%            ubar(i) ;
            switch i
               case 2
                Y = Zu(:,i);
                Yr = Zu(:,i+1);
                Yl = Zu(:,N+1);
               case N+1
                Y = Zu(:,i);
                Yr = Zu(:,2);
                Yl = Zu(:,i-1);
               otherwise
                Y = Zu(:,i);
                Yr = Zu(:,i+1);
                Yl = Zu(:,i-1);
               end
            if (p==2)
               
                  uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2);
                  uhr1 = Y(1) + Y(2)*(h(i)/2);
        
                  uhl1 = Y(1)+Y(2)*(-h(i)/2);
                  uhl2 = Yl(1) + Yl(2)*(h(i-1)/2);
                  
            elseif(p==4)
                uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3;
                uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3;
                uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3;
                uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3;
            elseif(p==5)
                uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3+Yr(5)*(-h(i+1)/2)^4;
                uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3 + Y(5)*(h(i)/2)^4;
                uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3+Y(5)*(-h(i)/2)^4;
                uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3 +Yl(5)*(h(i-1)/2)^4; 
            elseif(p==6)
                uhr2 = Yr(1)+Yr(2)*(-h(i+1)/2)+Yr(3)*(-h(i+1)/2)^2+Yr(4)*(-h(i+1)/2)^3+Yr(5)*(-h(i+1)/2)^4 +Yr(6)*(-h(i+1)/2)^5;
                uhr1 = Y(1) + Y(2)*(h(i)/2)+Y(3)*(h(i)/2)^2 + Y(4)*(h(i)/2)^3 + Y(5)*(h(i)/2)^4 +Y(6)*(h(i)/2)^5;
                uhl1 = Y(1)+Y(2)*(-h(i)/2) + Y(3)*(-h(i)/2)^2 + Y(4)*(-h(i)/2)^3+Y(5)*(-h(i)/2)^4 +Y(6)*(-h(i)/2)^5;
                uhl2 = Yl(1) + Yl(2)*(h(i-1)/2) + Yl(3)*(h(i-1)/2)^2 + Yl(4)*(h(i-1)/2)^3 +Yl(5)*(h(i-1)/2)^4 + Yl(6)*(h(i-1)/2)^5; 
            else
                error('2')
            end
            f(i) = f(i) + (ur1*(uhr1)-ul2*(uhl2))/h(i);

        end
        
             phi = -(ur1^2-ul2^2)/(2*h(i))-f(i);%burgers
             
%              if(p==2)
%                 ave = y(1) ;
%              elseif(p==3 || p==4)
%                  x2bi = (1/h(i))*( ((h(i)/2))^3/3-((-h(i)/2))^3/3);
%                 ave = y(1)+y(3)*x2bi;
%              elseif(p==5 || p==6)
%                 x2bi = (1/h(i))*( ((h(i)/2))^3/3-((-h(i)/2))^3/3);
%                 x4bi = (1/h(i))*( ((h(i)/2))^5/5-((-h(i)/2))^5/5);
%                 ave = y(1)+y(3)*x2bi+y(5)*x4bi;
%              end
%              if(ave<0)
%              phi = -(ur2^2-ul1^2)/(2*h(i))-f(i);%burgers
%              end
             
if( (ur1+ur2)<0 || (ul1+ul2)<0)
            fr = ur1;
            fl = ul2;
            if( (ur1+ur2)<0)
                fr = ur2;
            end
            if( (ul1+ul2)<0)
                fl = ul1;
            end
              phi = -(fr^2-fl^2)/(2*h(i))-f(i);%burgers
            
end
              
    else%residual eval
        sp = gsp(i);
        ut = fnval(fnder(sp),time);
    
        phi = -ut-(ur1^2-ul2^2)/(2*h(i))-f(i); % -(ur1*UU(i+1,j)-ul2*UU(i,j))/h(i);%burgers
        
      
        
    end
end


end

