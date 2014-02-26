function [ R ,uxx,Z] = computeres(u,x,k,h,N,f,r )
%COMPUTERES Summary of this function goes here
%   Detailed explanation goes here

switch r
case 4
    
    
    [err,Z]=unstructuredrecon3(u,x,h,N,NaN,NaN);
R = zeros(N+2,1);
uxx = zeros(N+2,1);
for i = 2:N+1
  %%%uxx(i) = 2*Z(3,i);
  
%   switch i
%     case 2
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,N+1);
%     case N+1
%         y = Z(:,i);
%         yr = Z(:,2);
%         yl = Z(:,i-1);
%     otherwise
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,i-1); 
% end    
%  upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2;
%         upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2;
%        
%         upr = (upr1+upr2)/2 ;%+(.2/((h(i)+h(i+1))/2))*(ur-ul);
% 
%         upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2;
%         upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2;
% 
%       
%         upl = (upl1+upl2)/2;% +(.2/((h(i)+h(i+1))/2))*(ur-ul);
%   uxx(i) = (upr-upl)/h(i);
%   
  
  
[ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r);%%%
 %%% R(i) = phi;%%%
  
%   R(i) = uxx(i) - f(i);
end



case 5
    [err,Z]=unstructuredrecon4(u,x,h,N,NaN,NaN);
R = zeros(N+2,1);
uxx = zeros(N+2,1);

for i = 2:N+1
%     hi = h(i);
%     xi = x(i);
%    x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);
%    
% %     uxx(i) = 2*Z(3,i)+12*Z(5,i)*x2bi;
% switch i
%     case 2
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,N+1);
%     case N+1
%         y = Z(:,i);
%         yr = Z(:,2);
%         yl = Z(:,i-1);
%     otherwise
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,i-1); 
% end    
% upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3;
% upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3;
% upr = (upr1+upr2)/2;
% upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3;
% upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3;
% upl = (upl1+upl2)/2;
%   uxx(i) = (upr-upl)/h(i);
  

   
%   R(i) = uxx(i) - f(i);

[ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r);%%%
  %%%R(i) = phi;%%%

end
    

case 6
    [err,Z]=unstructuredrecon5(u,x,h,N,NaN,NaN);
R = zeros(N+2,1);
uxx = zeros(N+2,1);

 for i = 2:N+1
%     hi = h(i);
%     xi = x(i);
%    x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);
%    
%     uxx(i) = 2*Z(3,i)+12*Z(5,i)*x2bi;

%  
%   switch i
%     case 2
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,N+1);
%     case N+1
%         y = Z(:,i);
%         yr = Z(:,2);
%         yl = Z(:,i-1);
%     otherwise
%         y = Z(:,i);
%         yr = Z(:,i+1);
%         yl = Z(:,i-1); 
% end    
%  upr1 = y(2)+2*y(3)*h(i)/2+3*y(4)*(h(i)/2)^2  + 4*y(5)*(h(i)/2)^3+5*y(6)*(h(i)/2)^4;
% upr2 = yr(2)+2*yr(3)*-h(i+1)/2+3*yr(4)*(-h(i+1)/2)^2   + 4*yr(5)*(-h(i+1)/2)^3 + 5*yr(6)*(-h(i+1)/2)^4;
% upr = (upr1+upr2)/2;
% upl1 = y(2)+2*y(3)*-h(i)/2+3*y(4)*(-h(i)/2)^2   + 4*y(5)*(-h(i)/2)^3+ 5*y(6)*(-h(i)/2)^4;
% upl2 = yl(2)+2*yl(3)*h(i-1)/2+3*yl(4)*(h(i-1)/2)^2   + 4*yl(5)*(h(i-1)/2)^3 + 5*yl(6)*(h(i-1)/2)^4;
% upl = (upl1+upl2)/2;
%   uxx(i) = (upr-upl)/h(i);

  

%   R(i) = uxx(i) - f(i);
[ur,ul,R(i)] = reconflux(u,Z,f,k,h,i,N,r);%%%
 %%% R(i) = phi;%%%

end
    


    otherwise 
               assert(0==1)
end
% [error,Z]=unstructuredrecon3(u,x,h,N,1,exp(-1));
% R = zeros(N+2,1);
% uxx = zeros(N+2,1);
% 
% for i = 2:N+1
%     hi = h(i);
%     xi = x(i);
%    %x2bi = (1/hi)*( ((xi+hi/2)-xi)^3/3-((xi-hi/2)-xi)^3/3);
%    
%     %uxx(i) = 2*Z(3,i)+12*Z(5,i)*x2bi;
%    
%   uxx(i) = 2*Z(3,i);
% %   uxx(i) = 2*Z(3,i);%+6*Z(4,i)*(x-
%    R(i) = uxx(i) - f(i);
% end
% uxx(N/2)
end

