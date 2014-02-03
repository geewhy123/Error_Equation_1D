function [ y ] = recon5( xi,hi,ui,x1,h1,u1, x2 ,h2,u2, x3,h3,u3,x4,h4,u4 ,x5,h5,u5,x6,h6,u6,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here




wi1 = 1;%1/abs(x1-xi);
wi2 = 1;%1/abs(x2-xi);
wi3 = 1;%1/abs(x3-xi);
wi4 = 1;%1/abs(x4-xi);
wi5 = 1;
wi6 = 1;


%xbi = (1/hi)*int(z-xi,z,xi-hi/2,xi+hi/2);
xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2);



%x2bi = (1/hi)*int((z-xi)^2,z,xi-hi/2,xi+hi/2);
x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);



%x3bi = (1/hi)*int((z-xi)^3,z,xi-hi/2,xi+hi/2);
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
% % A = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-x3bi); 
% %             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-x3bi); 
% %             wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) wi3*(x3b3+3*(x3-xi)*x2b3+3*(x3-xi)^2*xb3+(x3-xi)^3-x3bi);
% %             wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi) wi4*(x3b4+3*(x4-xi)*x2b4+3*(x4-xi)^2*xb4+(x4-xi)^3-x3bi)]);
 b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ; wi5*(ub5-ubi); wi6*(ub6-ubi)];

global AD
%%y(2:4) = A\b;%(A'*A)\(A'*b)

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

