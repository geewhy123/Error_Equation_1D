function [ upr,upl] = reconflux( u,Z,f,k,h,i,N,p)
%RECONFLUX Summary of this function goes here
%   Detailed explanation goes here

    y = Z(:,i);
    yr = Z(:,i+1);
    yl = Z(:,i-1);


switch p
    case 2
       upr1 = y(2);
upr2 = yr(2);
ur = yr(1)+yr(2)*(-h(i+1)/2);
ul = y(1) + y(2)*(h(i)/2);
upr = (upr1+upr2)/2 + (.2/h(i))*(ur(1)-ul(1));
upl1 = y(2);
upl2 = yl(2);

ur = y(1)+y(2)*(-h(i)/2);
ul = yl(1) + yl(2)*(h(i-1)/2);
upl = (upl1+upl2)/2 + (.2/h(i))*(ur(1)-ul(1));
        
        if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end



    case 3
        upr1 = y(2)+2*y(3)*h(i)/2; % u_i+1/2 using recon in i 
upr2 = yr(2)+2*yr(3)*-h(i+1)/2;
upr = (upr1+upr2)/2;
upl1 = y(2)+2*y(3)*-h(i)/2;
upl2 = yl(2)+2*yl(3)*h(i-1)/2;
upl = (upl1+upl2)/2;

        
   if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end


    case 4
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
        upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
        ur = yr(1)+yr(2)*(-h(i+1)/2)+yr(3)*(-h(i+1)/2)^2+yr(4)*(-h(i+1)/2)^3;
        ul = y(1) + y(2)*(h(i)/2)+y(3)*(h(i)/2)^2 + y(4)*(h(i)/2)^3;
        upr = (upr1+upr2)/2 ;

        upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
        upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;

        ur = y(1)+y(2)*(-h(i)/2) + y(3)*(-h(i)/2)^2 + y(4)*(-h(i)/2)^3;
        ul = yl(1) + yl(2)*(h(i-1)/2) + yl(3)*(h(i-1)/2)^2 + yl(4)*(h(i-1)/2)^3;
        upl = (upl1+upl2)/2 ;
        
if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end

    case 5
        upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3;
upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
upr = (upr1+upr2)/2;
upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
upl = (upl1+upl2)/2;
if i==2
    upl = upl1;
end
if i==N+1
   upr = upr1; 
end

        
        
    otherwise
        error('1')
end




end

