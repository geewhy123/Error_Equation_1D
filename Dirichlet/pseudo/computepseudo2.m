function [ AD ] = computepseudo2(N,x,h)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here

AD = zeros(2,4,N+2);
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
    h3 = h(i+2);%h(N+1);
    h4 = h(i+3);%h(i+2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i+2);%-1+x(N+1);
    x4 = x(i+3);%x(i+2);
    xi = x(i);
    end
    if i==N
    h1 = h(i-1);
    h2 = h(i+1);
    h3 = h(i-2);
    h4 = h(i-3);%h(2);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i+1);
    x3 = x(i-2);
    x4 = x(i-3);%1+x(2);
    xi = x(i); 
    end
    if i==2
         h1 = h(i+1);
    h2 = h(i+2);
    h3 = h(i+3);
    h4 = h(i+4);
    hi = h(i);
    x1 = x(i+1);
    x2 = x(i+2);
    x3 = x(i+3);
    x4 = x(i+4);
    xi = x(i);
    end
    if i==N+1
    h1 = h(i-1);
    h2 = h(i-2);
    h3 = h(i-3);
    h4 = h(i-4);
    hi = h(i);
    x1 = x(i-1);
    x2 = x(i-2);
    x3 = x(i-3);
    x4 = x(i-4);
    xi = x(i); 
    end
wi1 = 1;%/abs(x1-xi);
wi2 = 1;%/abs(x2-xi);
wi3 = 1;%/abs(x3-xi);
wi4 = 1;%/abs(x4-xi);
    
  
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



AA = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) ; 
             wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) ; 
             wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) ;
             wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi)  ]);
 AD(:,:,i)=pinv(AA);

 
%  Q = [-1 1 -2 2; 2.9412 3.9412 11.7647 11.7647];
%  if( (i>1)&&(i<N+2)&& (norm(AD(:,:,i)-Q) > 1e-2))
%      AD(:,:,i)
%     error('1') 
%  end
%  
%  if(i==4)
%      AA
%      error('1')
%  end
AA
end




end

