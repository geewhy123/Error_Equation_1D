function [ upr,upl,phi] = reconflux( u,Z,f,k,h,i,N,p)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here

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



ur = yr(1)+yr(2)*(-h(i+1)/2);
ul = y(1) + y(2)*(h(i)/2);
upr = (upr1+upr2)/2 + (.2/((h(i)+h(i+1))/2))*(ur-ul);
upl1 = y(2);
upl2 = yl(2);

UR = ur;

ur = y(1)+y(2)*(-h(i)/2);
ul = yl(1) + yl(2)*(h(i-1)/2);
upl = (upl1+upl2)/2 + (.2/((h(i)+h(i-1))/2))*(ur-ul);

UL = ur;


%        phi= (upr-upl)/h(i)-f(i);
 %       phi= (UR-UL)/h(i)-f(i);


    case 3
        upr1 = y(2)+2*y(3)*h(i)/2; % u_i+1/2 using recon in i 
upr2 = yr(2)+2*yr(3)*(-h(i+1)/2);
ur = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2;
ul = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2;

      
UR = ur;

upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur-ul);
upl1 = y(2)+2*y(3)*(-h(i)/2);
upl2 = yl(2)+2*yl(3)*h(i-1)/2;
 ur = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 ;
        ul = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 ;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ur-ul);

        
UL = ur;
        
   phi= (upr-upl)/h(i)-f(i);

    case 4
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
        upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
        ur = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3;
        ul = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3;
        upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur-ul);

                
UR = ur;

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;

        ur = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
        ul = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
        upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ur-ul);
     
        
UL = ur;
        
      phi= (upr-upl)/h(i)-f(i);

    case 5
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3;
upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
upr = (upr1+upr2)/2;

      ur = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4;
UR = ur;

upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
upl = (upl1+upl2)/2;

    ur = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+y(5)*(-h(i)/2)^4;
UL = ur;
   phi= (upr-upl)/h(i)-f(i);

    case 6
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3+5*y(6)*(h(i)/2)^4;
upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3 + 5*yr(6)*(-h(i+1)/2)^4;
upr = (upr1+upr2)/2;


UR = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3+yr(5)*(-h(i+1)/2)^4;


upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4;
upl = (upl1+upl2)/2;


UL = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3+ y(5)*(-h(i)/2)^4;
           phi= (upr-upl)/h(i)-f(i);
        
    otherwise
        error('1')
end



%        phi= (upr-upl)/h(i)-f(i);
        phi= (UR-UL)/h(i)-f(i);


end

