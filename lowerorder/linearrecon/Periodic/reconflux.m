function [ upr,upl,phi] = reconflux( u,Z,f,k,h,i,N,p,phys,uder,j,time,gsp,Rsp)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here

global TEND
if(abs(time-round(time))<1e-10)
   time = round(time); 
elseif(abs(time-TEND) < 1e-10)
   time = TEND;
end

global UU

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
        if(~isnan(time))
     
            %time
            sp = Rsp(i);
            f(i) = -1*getRes(time,k,i,sp);
        end
        phi= (ur2-ul1)/h(i)-f(i);
        
    else
%         if (i==5)
 %        uder(i,end)
%         error('1')
%         end


%%%phi= -uder(i,j)+(ur2-ul1)/h(i)-f(i);

sp = gsp(i);
ut = fnval(fnder(sp),time);


phi= -ut+(ur2-ul1)/h(i)-f(i);

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
    
    if((nargin < 12) || (isnan(time))|| isnan(j))
        phi = -(ur1^2-ul2^2)/(2*h(i))-f(i);%burgers
    else
        sp = gsp(i);
        ut = fnval(fnder(sp),time);
        
        phi = -ut-(ur1^2-ul2^2)/(2*h(i))-f(i) % -(ur1*UU(i+1,j)-ul2*UU(i,j))/h(i);%burgers
    end
end


end

