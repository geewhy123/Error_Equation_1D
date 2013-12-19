function [ AD ] = computepseudoshort(N,x,h)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here

AD = zeros(1,2,N+2);
for i = 2:N+1
    if i > 2 && i < N+1
    h1 = h(i-1);
    h2 = h(i+1);

    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);

    xi = x(i);
    end
%     if i==3
%          h1 = h(i-1);
%     h2 = h(i+1);
% 
%     hi = h(i);
%     x1 = x(i-1);
%     x2 = x(i+1);
% 
%     xi = x(i);
%     end
%     if i==N
%     h1 = h(i+1);
%     h2 = h(i-1);
%  
%     hi = h(i);
%     x1 = x(i+1);
%     x2 = x(i-1);
%  
%     xi = x(i); 
%     end
    if i==2
         h1 = h(i+1);
    h2 = h(i+2);

    hi = h(i);
    x1 = x(i+1);
    x2 = x(i+2);

    xi = x(i);
    end
    if i==N+1
    h1 = h(i-1);
    h2 = h(i-2);

    hi = h(i);
    x1 = x(i-1);
    x2 = x(i-2);

    xi = x(i); 
    end
wi1 = 1/abs(x1-xi);
wi2 = 1/abs(x2-xi);



xb1 = (1/h1)*( ((x1+h1/2)-x1)^2/2 -((x1-h1/2)-x1)^2/2 );
xb2 = (1/h2)*( ((x2+h2/2)-x2)^2/2 -((x2-h2/2)-x2)^2/2 );

xbi = (1/hi)*( ((xi+hi/2)-xi)^2/2 -((xi-hi/2)-xi)^2/2 );



AA = double([wi1*(xb1-xbi+x1-xi) ; 
            wi2*(xb2-xbi+x2-xi) ;  ]);
 AD(:,:,i)=pinv(AA);


end




end

