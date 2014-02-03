function [ AD ] = computepseudo1(N,x,h)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here

AD = zeros(1,4,N+2);
for i = 2:N+1
    if ((i > 3) && (i < N))
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
    

xb1 = (1/h1)*( ((x1+h1/2)-x1)^2/2 -((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*( ((x2+h2/2)-x2)^2/2 -((x2-h2/2)-x2)^2/2 );
xb3 = (1/h3)*( ((x3+h3/2)-x3)^2/2 -((x3-h3/2)-x3)^2/2 );
xb4 = (1/h4)*( ((x4+h4/2)-x4)^2/2 -((x4-h4/2)-x4)^2/2 );
xbi = (1/hi)*( ((xi+hi/2)-xi)^2/2 -((xi-hi/2)-xi)^2/2 );



AA = double([wi1*(xb1-xbi+x1-xi) ; 
            wi2*(xb2-xbi+x2-xi) ; 
            wi3*(xb3-xbi+x3-xi)  ;
            wi4*(xb4-xbi+x4-xi)  ]);
 AD(:,:,i)=pinv(AA);

 


end




end

