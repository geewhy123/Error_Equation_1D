function [ y ] = recon2( xi,hi,ui,x1,h1,u1, x2 ,h2,u2, x3,h3,u3,x4,h4,u4 ,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here




wi1 = 1;%1/abs(x1-xi);
wi2 = 1;%1/abs(x2-xi);
wi3 = 1;%1/abs(x3-xi);
wi4 = 1;%1/abs(x4-xi);


xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2);

x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);

ub1 = u1;
ub2 = u2;
ub3 = u3;
ub4 = u4;
ubi = ui;
 b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];

global AD

y(2:3) = AD(:,:,i)*b;

%y-y2;

y(1) = ubi-xbi*y(2)-x2bi*y(3);%ubi-xbi*y(2)
%q = y(1)-ubi

end

