function [ y ] = recon1( xi,hi,ui,x1,h1,u1, x2 ,h2,u2, x3,h3,u3,x4,h4,u4 ,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here




wi1 = 1/abs(x1-xi)^0;
wi2 = 1/abs(x2-xi)^0;
wi3 = 1/abs(x3-xi)^0;
wi4 = 1/abs(x4-xi)^0;



% xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2);

ub1 = u1;
ub2 = u2;
ub3 = u3;
ub4 = u4;
ubi = ui;

b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];

global AD


y(2) = AD(:,:,i)*b;


y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
%q = y(1)-ubi


end

