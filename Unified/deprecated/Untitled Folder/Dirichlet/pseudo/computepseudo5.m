function [ AD ] = computepseudo5(N,x,h)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here

AD = zeros(5,6,N+2);
for i = 2:N+1
    if i > 4 && i < N-1
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i+2);
    h5 = h(i-3);
    h6 = h(i+3);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i+2);
    x5 = x(i-3);
    x6 = x(i+3);
    xi = x(i);
    end
    if i==4
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i+2);
    h5 = h(i+3);
    h6 = h(i+4);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i+2);
    x5 = x(i+3);%-1+x(N+1);
    x6 = x(i+4);%x(i+3);
    xi = x(i);
    end
    if i==N-1
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i+2);
    h5 = h(i-3);
    h6 = h(i-4);%h(2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i+2);
    x5 = x(i-3);
    x6 = x(i-4);%1+x(2);
    xi = x(i);
    end
     if i==3
         h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i+2);%h(N+1);
    h4 = h(i+3);%h(i+2);
    h5 = h(i+4);%h(N);
    h6 = h(i+5);%h(i+3);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i+2);%-1+x(N+1);
    x4 = x(i+3);%x(i+2);
    x5 = x(i+4);%-1+x(N);
    x6 = x(i+5);%x(i+3);
    xi = x(i);
    end
    if i==N
   h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i-3);%h(2);
    h5 = h(i-4);%h(i-3);
    h6 = h(i-5);%h(3);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i-3);%1+x(2);
    x5 = x(i-4);%x(i-3);
    x6 = x(i-5);%1+x(3);
    xi = x(i);
    end
    if i==2
%     h1 = h(N+1);
%     h2 = h(i+1);
%     h3 = h(N);
%     h4 = h(i+2);
%     h5 = h(N-1);
%     h6 = h(i+3);
    h1 = h(i+1);
    h2 = h(i+2);
    h3 = h(i+3);
    h4 = h(i+4);
    h5 = h(i+5);
    h6 = h(i+6);
    hi = h(i);
%     x1 = -1+x(N+1);
%     x2 = x(i+1);
%     x3 = -1+x(N);
%     x4 = x(i+2);
%     x5 = -1+x(N-1);
%     x6 = x(i+3);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);
    x4 = x(i+4);
    x5 = x(i+5);
    x6 = x(i+6);

    xi = x(i);
    end
    if i==N+1
     h1 = h(i-1);
    h2 = h(i-2);
    h3 = h(i-3);
    h4 = h(i-4);
    h5 = h(i-5);
    h6 = h(i-6);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i-2);
    x3 = x(i-3);
    x4 = x(i-4);
    x5 = x(i-5);
    x6 = x(i-6);
    xi = x(i); 
    end
    wi1 = 1;
    wi2 = 1;
    wi3 = 1;
    wi4 = 1;
    wi5 = 1;
    wi6 = 1;
    
xb1 = (1/h1)*( ((x1+h1/2)-x1)^2/2 -((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*( ((x2+h2/2)-x2)^2/2 -((x2-h2/2)-x2)^2/2 );
xb3 = (1/h3)*( ((x3+h3/2)-x3)^2/2 -((x3-h3/2)-x3)^2/2 );
xb4 = (1/h4)*( ((x4+h4/2)-x4)^2/2 -((x4-h4/2)-x4)^2/2 );
xb5 = (1/h5)*( ((x5+h5/2)-x5)^2/2 -((x5-h5/2)-x5)^2/2 );
xb6 = (1/h6)*( ((x6+h6/2)-x6)^2/2 -((x6-h6/2)-x6)^2/2 );
xbi = (1/hi)*( ((xi+hi/2)-xi)^2/2 -((xi-hi/2)-xi)^2/2 );

x2b1 = (1/h1)*( ((x1+h1/2)-x1)^3/3 -((x1-h1/2)-x1)^3/3 );
x2b2 = (1/h2)*( ((x2+h2/2)-x2)^3/3 -((x2-h2/2)-x2)^3/3 );
x2b3 = (1/h3)*( ((x3+h3/2)-x3)^3/3 -((x3-h3/2)-x3)^3/3 );
x2b4 = (1/h4)*( ((x4+h4/2)-x4)^3/3 -((x4-h4/2)-x4)^3/3 );
x2b5 = (1/h5)*( ((x5+h5/2)-x5)^3/3 -((x5-h5/2)-x5)^3/3 );
x2b6 = (1/h6)*( ((x6+h6/2)-x6)^3/3 -((x6-h6/2)-x6)^3/3 );
x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3 -((xi-hi/2)-xi)^3/3 );

x3b1 = (1/h1)*( ((x1+h1/2)-x1)^4/4 -((x1-h1/2)-x1)^4/4 );
x3b2 = (1/h2)*( ((x2+h2/2)-x2)^4/4 -((x2-h2/2)-x2)^4/4 );
x3b3 = (1/h3)*( ((x3+h3/2)-x3)^4/4 -((x3-h3/2)-x3)^4/4 );
x3b4 = (1/h4)*( ((x4+h4/2)-x4)^4/4 -((x4-h4/2)-x4)^4/4 );
x3b5 = (1/h5)*( ((x5+h5/2)-x5)^4/4 -((x5-h5/2)-x5)^4/4 );
x3b6 = (1/h6)*( ((x6+h6/2)-x6)^4/4 -((x6-h6/2)-x6)^4/4 );
x3bi = (1/hi)*( ((xi+hi/2)-xi)^4/4 -((xi-hi/2)-xi)^4/4 );

x4b1 = (1/h1)*( ((x1+h1/2)-x1)^5/5 -((x1-h1/2)-x1)^5/5 );
x4b2 = (1/h2)*( ((x2+h2/2)-x2)^5/5 -((x2-h2/2)-x2)^5/5 );
x4b3 = (1/h3)*( ((x3+h3/2)-x3)^5/5 -((x3-h3/2)-x3)^5/5 );
x4b4 = (1/h4)*( ((x4+h4/2)-x4)^5/5 -((x4-h4/2)-x4)^5/5 );
x4b5 = (1/h5)*( ((x5+h5/2)-x5)^5/5 -((x5-h5/2)-x5)^5/5 );
x4b6 = (1/h6)*( ((x6+h6/2)-x6)^5/5 -((x6-h6/2)-x6)^5/5 );
x4bi = (1/hi)*( ((xi+hi/2)-xi)^5/5 -((xi-hi/2)-xi)^5/5 );

x5b1 = (1/h1)*( ((x1+h1/2)-x1)^6/6 -((x1-h1/2)-x1)^6/6 );
x5b2 = (1/h2)*( ((x2+h2/2)-x2)^6/6 -((x2-h2/2)-x2)^6/6 );
x5b3 = (1/h3)*( ((x3+h3/2)-x3)^6/6 -((x3-h3/2)-x3)^6/6 );
x5b4 = (1/h4)*( ((x4+h4/2)-x4)^6/6 -((x4-h4/2)-x4)^6/6 );
x5b5 = (1/h5)*( ((x5+h5/2)-x5)^6/6 -((x5-h5/2)-x5)^6/6 );
x5b6 = (1/h6)*( ((x6+h6/2)-x6)^6/6 -((x6-h6/2)-x6)^6/6 );
x5bi = (1/hi)*( ((xi+hi/2)-xi)^6/6 -((xi-hi/2)-xi)^6/6 );


AA = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) wi1*(x3b1+3*(x1-xi)*x2b1+3*(x1-xi)^2*xb1+(x1-xi)^3-x3bi)   wi1*(x4b1+4*(x1-xi)*x3b1+6*(x1-xi)^2*x2b1+4*(x1-xi)^3*xb1+(x1-xi)^4-x4bi) wi1*(x5b1+5*(x1-xi)*x4b1+10*(x1-xi)^2*x3b1+10*(x1-xi)^3*x2b1+5*(x1-xi)^4*xb1+(x1-xi)^5-x5bi); 
             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) wi2*(x3b2+3*(x2-xi)*x2b2+3*(x2-xi)^2*xb2+(x2-xi)^3-x3bi)   wi2*(x4b2+4*(x2-xi)*x3b2+6*(x2-xi)^2*x2b2+4*(x2-xi)^3*xb2+(x2-xi)^4-x4bi) wi2*(x5b2+5*(x2-xi)*x4b2+10*(x2-xi)^2*x3b2+10*(x2-xi)^3*x2b2+5*(x2-xi)^4*xb2+(x2-xi)^5-x5bi); 
             wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) wi3*(x3b3+3*(x3-xi)*x2b3+3*(x3-xi)^2*xb3+(x3-xi)^3-x3bi)   wi3*(x4b3+4*(x3-xi)*x3b3+6*(x3-xi)^2*x2b3+4*(x3-xi)^3*xb3+(x3-xi)^4-x4bi) wi3*(x5b3+5*(x3-xi)*x4b3+10*(x3-xi)^2*x3b3+10*(x3-xi)^3*x2b3+5*(x3-xi)^4*xb3+(x3-xi)^5-x5bi);
             wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi) wi4*(x3b4+3*(x4-xi)*x2b4+3*(x4-xi)^2*xb4+(x4-xi)^3-x3bi)   wi4*(x4b4+4*(x4-xi)*x3b4+6*(x4-xi)^2*x2b4+4*(x4-xi)^3*xb4+(x4-xi)^4-x4bi) wi4*(x5b4+5*(x4-xi)*x4b4+10*(x4-xi)^2*x3b4+10*(x4-xi)^3*x2b4+5*(x4-xi)^4*xb4+(x4-xi)^5-x5bi); 
             wi5*(xb5-xbi+x5-xi) wi5*(x2b5+2*(x5-xi)*xb5+(x5-xi)^2-x2bi) wi5*(x3b5+3*(x5-xi)*x2b5+3*(x5-xi)^2*xb5+(x5-xi)^3-x3bi)   wi5*(x4b5+4*(x5-xi)*x3b5+6*(x5-xi)^2*x2b5+4*(x5-xi)^3*xb5+(x5-xi)^4-x4bi) wi5*(x5b5+5*(x5-xi)*x4b5+10*(x5-xi)^2*x3b5+10*(x5-xi)^3*x2b5+5*(x5-xi)^4*xb5+(x5-xi)^5-x5bi);
             wi6*(xb6-xbi+x6-xi) wi6*(x2b6+2*(x6-xi)*xb6+(x6-xi)^2-x2bi) wi6*(x3b6+3*(x6-xi)*x2b6+3*(x6-xi)^2*xb6+(x6-xi)^3-x3bi)   wi6*(x4b6+4*(x6-xi)*x3b6+6*(x6-xi)^2*x2b6+4*(x6-xi)^3*xb6+(x6-xi)^4-x4bi) wi6*(x5b6+5*(x6-xi)*x4b6+10*(x6-xi)^2*x3b6+10*(x6-xi)^3*x2b6+5*(x6-xi)^4*xb6+(x6-xi)^5-x5bi);]);
 AD(:,:,i)=pinv(AA);


end




end

