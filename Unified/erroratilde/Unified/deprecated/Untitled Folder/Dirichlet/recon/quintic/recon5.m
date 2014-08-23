function [ y ] = recon5( xi,hi,ui,x1,h1,u1, x2 ,h2,u2, x3,h3,u3,x4,h4,u4 ,x5,h5,u5,x6,h6,u6,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here




wi1 = 1;%1/abs(x1-xi);
wi2 = 1;%1/abs(x2-xi);
wi3 = 1;%1/abs(x3-xi);
wi4 = 1;%1/abs(x4-xi);
wi5 = 1;
wi6 = 1;


xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2);
x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);
x3bi = (1/hi)*( ((xi+hi/2)-xi)^4/4-((xi-hi/2)-xi)^4/4);
x4bi = (1/hi)*( ((xi+hi/2)-xi)^5/5-((xi-hi/2)-xi)^5/5);
x5bi = (1/hi)*( ((xi+hi/2)-xi)^6/6-((xi-hi/2)-xi)^6/6);



ub1 = u1;
ub2 = u2;
ub3 = u3;
ub4 = u4;
ub5 = u5;
ub6 = u6;
ubi = ui;
 b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ; wi5*(ub5-ubi); wi6*(ub6-ubi)];

global AD


y(2:6) = AD(:,:,i)*b;

%y-y2;

y(1) = ubi-xbi*y(2)-x2bi*y(3)-x3bi*y(4)-x4bi*y(5)-x5bi*y(6);%ubi-xbi*y(2)
%q = y(1)-ubi

%if xi==0.1
%     wi1
%     x2b1
%     x1
%     xi
    %y
%     x2bi
%     y
%   x2b2 
%end

end

