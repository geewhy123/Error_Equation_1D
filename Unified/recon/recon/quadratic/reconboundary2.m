function [ y ] = reconboundary2( xi,hi,ui,x1,h1,u1, x2 ,h2,u2, x3,h3,u3,x4,h4,u4 ,uL,str,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here

%global AD


wi1 = 1;%/abs(x1-xi)^w;
wi2 = 1;%/abs(x2-xi)^w;
wi3 = 1;%/abs(x3-xi)^w;
wi4 = 1;%/abs(x4-xi)^w;

xb1 = (1/h1)*(((x1+h1/2)-x1)^2/2-((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*(((x2+h2/2)-x2)^2/2-((x2-h2/2)-x2)^2/2 );
xb3 = (1/h3)*(((x3+h3/2)-x3)^2/2-((x3-h3/2)-x3)^2/2 );
xb4 = (1/h4)*(((x4+h4/2)-x4)^2/2-((x4-h4/2)-x4)^2/2 );
xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2 );

x2b1 = (1/h1)*(((x1+h1/2)-x1)^3/3-((x1-h1/2)-x1)^3/3 );
x2b2 = (1/h2)*(((x2+h2/2)-x2)^3/3-((x2-h2/2)-x2)^3/3 );
x2b3 = (1/h3)*(((x3+h3/2)-x3)^3/3-((x3-h3/2)-x3)^3/3 );
x2b4 = (1/h4)*(((x4+h4/2)-x4)^3/3-((x4-h4/2)-x4)^3/3 );
x2bi = (1/hi)*(((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3 );

x3b1 = (1/h1)*(((x1+h1/2)-x1)^4/4-((x1-h1/2)-x1)^4/4 );
x3b2 = (1/h2)*(((x2+h2/2)-x2)^4/4-((x2-h2/2)-x2)^4/4 );
x3b3 = (1/h3)*(((x3+h3/2)-x3)^4/4-((x3-h3/2)-x3)^4/4 );
x3b4 = (1/h4)*(((x4+h4/2)-x4)^4/4-((x4-h4/2)-x4)^4/4 );
x3bi = (1/hi)*(((xi+hi/2)-xi)^4/4-((xi-hi/2)-xi)^4/4 );

x4b1 = (1/h1)*(((x1+h1/2)-x1)^5/5-((x1-h1/2)-x1)^5/5 );
x4b2 = (1/h2)*(((x2+h2/2)-x2)^5/5-((x2-h2/2)-x2)^5/5 );
x4b3 = (1/h3)*(((x3+h3/2)-x3)^5/5-((x3-h3/2)-x3)^5/5 );
x4b4 = (1/h4)*(((x4+h4/2)-x4)^5/5-((x4-h4/2)-x4)^5/5 );
x4bi = (1/hi)*(((xi+hi/2)-xi)^5/5-((xi-hi/2)-xi)^5/5 );

ub1 = u1;
ub2 = u2;
ub3 = u3;
ub4 = u4;
ubi = ui;
A = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-x3bi)  wi1*(x4b1+4*(x1-xi)*x3b1+6*(x1-xi)^2*x2b1+4*(x1-xi)^3*xb1+(x1-xi)^4-x4bi); 
            wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-x3bi)  wi2*(x4b2+4*(x2-xi)*x3b2+6*(x2-xi)^2*x2b2+4*(x2-xi)^3*xb2+(x2-xi)^4-x4bi); 
            wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) wi3*(x3b3+3*(x3-xi)*x2b3+3*(x3-xi)^2*xb3+(x3-xi)^3-x3bi)  wi3*(x4b3+4*(x3-xi)*x3b3+6*(x3-xi)^2*x2b3+4*(x3-xi)^3*xb3+(x3-xi)^4-x4bi);
            wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi) wi4*(x3b4+3*(x4-xi)*x2b4+3*(x4-xi)^2*xb4+(x4-xi)^3-x3bi)  wi4*(x4b4+4*(x4-xi)*x3b4+6*(x4-xi)^2*x2b4+4*(x4-xi)^3*xb4+(x4-xi)^4-x4bi)]);
b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];


%A = AD(:,:,1);

%uL = 1;

if strcmp(str,'right')
   hi = -hi; 
end

A = (A(:,2:2)-([(wi1*(xb1+x1-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2) ; 
                (wi2*(xb2+x2-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2)  ;
                (wi3*(xb3+x3-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2)  ;
                (wi4*(xb4+x4-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2) ]));
   
                 
b = (b-[(wi1*(xb1+x1-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ; 
       (wi2*(xb2+x2-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ;
       (wi3*(xb3+x3-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ;
       (wi4*(xb4+x4-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ;]);

A
b
error('1')

y(3:3) = (A'*A)\(A'*b);

P = [ 1 -hi/2; 1 xbi];
q = [uL; ubi]-[ (-hi/2)^2*y(3) ; x2bi*y(3)];
y(1:2) = P\q;
%y = double(y);

%y(1) = ubi-xbi*y(2)-x2bi*y(3)-x3bi*y(4);%ubi-xbi*y(2)

%y(1) = uL-(-h1/2)*y(2)-(-h1/2)^2*y(3)-(-h1/2)^3*y(4)
double(y(1)+xbi*y(2)+x2bi*y(3));
double(y(1)+y(2)*(-hi/2)+y(3)*(-hi/2)^2);

%q = y(1)-ubi


end

