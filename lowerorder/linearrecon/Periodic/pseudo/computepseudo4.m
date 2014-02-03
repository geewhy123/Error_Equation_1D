function [ AD ] = computepseudo4(N,x,h)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here

AD = zeros(4,4,N+2);
for i = 2:N+1
    if i > 3 && i < N
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i+2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i+2);
    xi = x(i);
    end
     if i==3
         h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(N+1);
    h4 = h(i+2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = -1+x(N+1);
    x4 = x(i+2);
    xi = x(i);
    end
    if i==N
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = 1+x(2);
    xi = x(i); 
    end
    if i==2
         h1 = h(N+1);
    h2 = h(i+1);
    h3 = h(N);
    h4 = h(i+2);
    hi = h(i);
    x1 = -1+x(N+1);
    x2 = x(i+1);
    x3 = -1+x(N);
    x4 = x(i+2);
    xi = x(i);
    end
    if i==N+1
    h1 = h(i-1);
    h2 = h(2);
    h3 = h(i-2);
    h4 = h(3);
    hi = h(i);
    x1 = x(i-1);
    x2 = 1+x(2);
    x3 = x(i-2);
    x4 = 1+x(3);
    xi = x(i); 
    end
    wi1 = 1;
    wi2 = 1;
    wi3 = 1;
    wi4 = 1;
    
    syms z
xb1 = (1/h1)*( ((x1+h1/2)-x1)^2/2 -((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*( ((x2+h2/2)-x2)^2/2 -((x2-h2/2)-x2)^2/2 );
xb3 = (1/h3)*( ((x3+h3/2)-x3)^2/2 -((x3-h3/2)-x3)^2/2 );
xb4 = (1/h4)*( ((x4+h4/2)-x4)^2/2 -((x4-h4/2)-x4)^2/2 );
xbi = (1/hi)*( ((xi+hi/2)-xi)^2/2 -((xi-hi/2)-xi)^2/2 );

x2b1 = (1/h1)*( ((x1+h1/2)-x1)^3/3 -((x1-h1/2)-x1)^3/3 );
x2b2 = (1/h2)*( ((x2+h2/2)-x2)^3/3 -((x2-h2/2)-x2)^3/3 );
x2b3 = (1/h3)*( ((x3+h3/2)-x3)^3/3 -((x3-h3/2)-x3)^3/3 );
x2b4 = (1/h4)*( ((x4+h4/2)-x4)^3/3 -((x4-h4/2)-x4)^3/3 );
x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3 -((xi-hi/2)-xi)^3/3 );

x3b1 = (1/h1)*( ((x1+h1/2)-x1)^4/4 -((x1-h1/2)-x1)^4/4 );
x3b2 = (1/h2)*( ((x2+h2/2)-x2)^4/4 -((x2-h2/2)-x2)^4/4 );
x3b3 = (1/h3)*( ((x3+h3/2)-x3)^4/4 -((x3-h3/2)-x3)^4/4 );
x3b4 = (1/h4)*( ((x4+h4/2)-x4)^4/4 -((x4-h4/2)-x4)^4/4 );
x3bi = (1/hi)*( ((xi+hi/2)-xi)^4/4 -((xi-hi/2)-xi)^4/4 );

x4b1 = (1/h1)*( ((x1+h1/2)-x1)^5/5 -((x1-h1/2)-x1)^5/5 );
x4b2 = (1/h2)*( ((x2+h2/2)-x2)^5/5 -((x2-h2/2)-x2)^5/5 );
x4b3 = (1/h3)*( ((x3+h3/2)-x3)^5/5 -((x3-h3/2)-x3)^5/5 );
x4b4 = (1/h4)*( ((x4+h4/2)-x4)^5/5 -((x4-h4/2)-x4)^5/5 );
x4bi = (1/hi)*( ((xi+hi/2)-xi)^5/5 -((xi-hi/2)-xi)^5/5 );

AA = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-x3bi)  wi1*(x4b1+4*(x1-xi)*x3b1+6*(x1-xi)^2*x2b1+4*(x1-xi)^3*xb1+(x1-xi)^4-x4bi); 
             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-x3bi)   wi2*(x4b2+4*(x2-xi)*x3b2+6*(x2-xi)^2*x2b2+4*(x2-xi)^3*xb2+(x2-xi)^4-x4bi); 
             wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) wi3*(x3b3+3*(x3-xi)*x2b3+3*(x3-xi)^2*xb3+(x3-xi)^3-x3bi)   wi3*(x4b3+4*(x3-xi)*x3b3+6*(x3-xi)^2*x2b3+4*(x3-xi)^3*xb3+(x3-xi)^4-x4bi);
             wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi) wi4*(x3b4+3*(x4-xi)*x2b4+3*(x4-xi)^2*xb4+(x4-xi)^3-x3bi)   wi4*(x4b4+4*(x4-xi)*x3b4+6*(x4-xi)^2*x2b4+4*(x4-xi)^3*xb4+(x4-xi)^4-x4bi) ]);
 AD(:,:,i)=pinv(AA);

 if(i==4)
3*(x1-xi)*x2bi
(x1-xi)^3

AA
 end
 
end

error('1')



end

