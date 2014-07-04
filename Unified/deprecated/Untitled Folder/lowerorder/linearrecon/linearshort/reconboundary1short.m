function [ y ] = reconboundary1short( xi,hi,ui,x1,h1,u1, x2 ,h2,u2 ,uL,str,i)
%RECON3 Summary of this function goes here
%   Detailed explanation goes here

%global AD


wi1 = 1/abs(x1-xi);
wi2 = 1/abs(x2-xi);

% xi = x(i);
% x1 = x(i+1);
% x2 = x(i+2);
% x3 = x(i+3);
% x4 = x(i+4);
%syms z
% xb1 = (1/h1)*int(z-x1,z,x1-h1/2,x1+h1/2);
% xb2 = (1/h2)*int(z-x2,z,x2-h2/2,x2+h2/2);
% xb3 = (1/h3)*int(z-x3,z,x3-h3/2,x3+h3/2);
% xb4 = (1/h4)*int(z-x4,z,x4-h4/2,x4+h4/2);
% xbi = (1/hi)*int(z-xi,z,xi-hi/2,xi+hi/2);
% 
% x2b1 = (1/h1)*int((z-x1)^2,z,x1-h1/2,x1+h1/2);
% x2b2 = (1/h2)*int((z-x2)^2,z,x2-h2/2,x2+h2/2);
% x2b3 = (1/h3)*int((z-x3)^2,z,x3-h3/2,x3+h3/2);
% x2b4 = (1/h4)*int((z-x4)^2,z,x4-h4/2,x4+h4/2);
% x2bi = (1/hi)*int((z-xi)^2,z,xi-hi/2,xi+hi/2);
% 
% x3b1 = (1/h1)*int((z-x1)^3,z,x1-h1/2,x1+h1/2);
% x3b2 = (1/h2)*int((z-x2)^3,z,x2-h2/2,x2+h2/2);
% x3b3 = (1/h3)*int((z-x3)^3,z,x3-h3/2,x3+h3/2);
% x3b4 = (1/h4)*int((z-x4)^3,z,x4-h4/2,x4+h4/2);
% x3bi = (1/hi)*int((z-xi)^3,z,xi-hi/2,xi+hi/2);
xb1 = (1/h1)*(((x1+h1/2)-x1)^2/2-((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*(((x2+h2/2)-x2)^2/2-((x2-h2/2)-x2)^2/2 );

xbi = (1/hi)*(((xi+hi/2)-xi)^2/2-((xi-hi/2)-xi)^2/2 );

x2b1 = (1/h1)*(((x1+h1/2)-x1)^3/3-((x1-h1/2)-x1)^3/3 );
x2b2 = (1/h2)*(((x2+h2/2)-x2)^3/3-((x2-h2/2)-x2)^3/3 );

x2bi = (1/hi)*(((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3 );

% x3b1 = (1/h1)*(((x1+h1/2)-x1)^4/4-((x1-h1/2)-x1)^4/4 );
% x3b2 = (1/h2)*(((x2+h2/2)-x2)^4/4-((x2-h2/2)-x2)^4/4 );
% 
% x3bi = (1/hi)*(((xi+hi/2)-xi)^4/4-((xi-hi/2)-xi)^4/4 );
% 
% x4b1 = (1/h1)*(((x1+h1/2)-x1)^5/5-((x1-h1/2)-x1)^5/5 );
% x4b2 = (1/h2)*(((x2+h2/2)-x2)^5/5-((x2-h2/2)-x2)^5/5 );
% 
% x4bi = (1/hi)*(((xi+hi/2)-xi)^5/5-((xi-hi/2)-xi)^5/5 );
% 
% ub1 = u1;
% ub2 = u2;
% 
 ubi = ui;
% A = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-x3bi)  wi1*(x4b1+4*(x1-xi)*x3b1+6*(x1-xi)^2*x2b1+4*(x1-xi)^3*xb1+(x1-xi)^4-x4bi); 
%             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-x3bi)  wi2*(x4b2+4*(x2-xi)*x3b2+6*(x2-xi)^2*x2b2+4*(x2-xi)^3*xb2+(x2-xi)^4-x4bi); ]);
% b = [wi1*(ub1-ubi); wi2*(ub2-ubi);  ];


%A = AD(:,:,1);

%uL = 1;

if strcmp(str,'right')
   hi = -hi; 
end

% A = (A(:,2:2)-([(wi1*(xb1+x1-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2) ; 
%                      (wi2*(xb2+x2-xi-xbi)/(xbi-(-hi/2)))*(x2bi-(-hi/2)^2)  ;
%                  ]));
%    
%                  
% b = (b-[(wi1*(xb1+x1-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ; 
%        (wi2*(xb2+x2-xi-xbi)/(xbi-(-hi/2)))*(ubi-uL) ;
%   ;]);



%y(3:3) = (A'*A)\(A'*b);

P = [ 1 -hi/2; 1 xbi];
%q = [uL; ubi]-[ (-hi/2)^2*y(3) ; x2bi*y(3)];
q = [uL;ubi];
y(1:2) = P\q;
%y = double(y);

%y(1) = ubi-xbi*y(2)-x2bi*y(3)-x3bi*y(4);%ubi-xbi*y(2)

%y(1) = uL-(-h1/2)*y(2)-(-h1/2)^2*y(3)-(-h1/2)^3*y(4)
double(y(1)+xbi*y(2));
double(y(1)+y(2)*(-hi/2));


    
%q = y(1)-ubi


end

