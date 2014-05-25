function [ Z] = unstructuredrecon(obj,u,p,eqn)
%UNSTRUCTUREDRECON3 Summary of this function goes here
% p = obj.pOrder;
switch p
%     case 2
%    [Z] = unstructuredrecon1 (obj,u,p); 
%     case 3
%         [Z] = unstructuredrecon1 (obj,u,p); 
% %    [~,Z] = unstructuredrecon2 (u,x,h,N,u0,u1); 
%     case 4
%    [~,Z] = unstructuredrecon3 (u,x,h,N,u0,u1); 
%     case 5
%    [~,Z] = unstructuredrecon4 (u,x,h,N,u0,u1); 
        case 6
   [~,Z] = unstructuredrecon5 (u,x,h,N,u0,u1); 
    otherwise 
        [Z] = unstructuredreconp (obj,u,p,eqn); 
%         p
%         
%        assert(0==1)
end


end


function [Z] = unstructuredreconp(obj,u ,p,eqn)


%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
N = obj.nCells;
x = obj.cellCentroids;
h = obj.cellWidths;
moments = obj.moments;
Z = zeros(p,N+2);
% error = 0;


if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
    if(strcmp(eqn,'solution')==1)
       AD = obj.primalPI;
    elseif(strcmp(eqn,'error')==1)
%         obj.bcLeftType = 'D';
%         obj.bcLeftVal = 1;
%         obj.bcRightType = 'D';
%         obj.bcRightVal = 1;
       AD = obj.errorPI;
    elseif(strcmp(eqn,'residual')==1)
       AD = obj.resPI;
    else
        eqn
        error('1')
end
    
    
for i = 2:N+1
switch i
    case 2
cv1 = N+1;
cv2 = i+1;
cv3 = N;
cv4 = i+2;
    case 3
        cv1 = i-1;
cv2 = i+1;
cv3 = N+1;
cv4 = i+2;
    case N
        cv1 = i-1;
cv2 = i+1;
cv3 = i-2;
cv4 = 2;
    case N+1
        cv1 = i-1;
cv2 = 2;
cv3 = i-2;
cv4 = 3;
    otherwise
        cv1 =i-1;
cv2 = i+1;
cv3 = i-2;
cv4 = i+2;
 
end

wi1 = 1/abs(x(cv1)-x(i))^0;
wi2 = 1/abs(x(cv2)-x(i))^0;
wi3 = 1/abs(x(cv3)-x(i))^0;
wi4 = 1/abs(x(cv4)-x(i))^0;

ub1 = u(cv1);
ub2 = u(cv2);
ub3 = u(cv3);
ub4 = u(cv4);
ubi = u(i);


b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];

% AD = obj.primalPI;

Y(2:p) = AD(:,:,i)*b;


Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
%q = y(1)-ubi


for k = 1:p-1
   Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
end

 Z(:,i) = Y;
 
 
%  b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];


% y(1) = ubi-xbi*y(2)-x2bi*y(3);%ubi-xbi*y(2)
%q = y(1)-ubi

 
 
 
 
end

elseif(obj.bcLeftType=='D' && obj.bcRightType == 'D')
    
%     AA = obj.reconM;
% AA = zeros(4,p-1,N+2);
  if(strcmp(eqn,'solution')==1)
       AA = obj.primalRM;
       AD = obj.primalPI;
    elseif(strcmp(eqn,'error')==1)
%         obj.bcLeftType = 'D';
%         obj.bcLeftVal = 1;
%         obj.bcRightType = 'D';
%         obj.bcRightVal = 1;
       AA = obj.errorRM;
       AD = obj.errorPI;
    elseif(strcmp(eqn,'residual')==1)
       AA = obj.resRM;
       AD = obj.resPI;
    else
        eqn
        error('1')
end
    
    
i=2;
cv1 =i+1;
cv2 = i+2;
cv3 = i+3;
cv4 = i+4;
wi1 = 1/abs(x(cv1)-x(i))^0;
wi2 = 1/abs(x(cv2)-x(i))^0;
wi3 = 1/abs(x(cv3)-x(i))^0;
wi4 = 1/abs(x(cv4)-x(i))^0;

ub1 = u(cv1);
ub2 = u(cv2);
ub3 = u(cv3);
ub4 = u(cv4);
ubi = u(i);
uL = obj.bcLeftVal;



A = AA(:,:,i);
xbi = obj.moments(i,2);
x2bi = obj.moments(i,3);
if(p>2)
 b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
 
b = (b-[(A(1,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ; 
       (A(2,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
       (A(3,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
       (A(4,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;]);
    
%    
% C =           ([(A(1,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2) ; 
%                 (A(2,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2)  ;
%                 (A(3,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2)  ;
%                 (A(4,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2) ]);   
C = zeros(4,p-2);
for jj = 1:p-2
   C(1,jj) = (A(1,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1)); 
   C(2,jj) = (A(2,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
   C(3,jj) = (A(3,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
   C(4,jj) = (A(4,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
end

A = (A(:,2:p-1)-C);
        

y(3:p) = (A'*A)\(A'*b);
end
P = [ 1 -h(i)/2; 1 xbi];
q = [uL; ubi];
for j = 2:p-1
q = q-[ (-h(i)/2)^j*y(j+1) ; obj.moments(i,j+1)*y(j+1)];
end




y(1:2) = P\q;

Z(:,i) = y;


i = N+1;
xbi = obj.moments(i,2);
x2bi = obj.moments(i,3);
A = AA(:,:,i);
cv1 =i-1;
cv2 = i-2;
cv3 = i-3;
cv4 = i-4;
wi1 = 1/abs(x(cv1)-x(i))^0;
wi2 = 1/abs(x(cv2)-x(i))^0;
wi3 = 1/abs(x(cv3)-x(i))^0;
wi4 = 1/abs(x(cv4)-x(i))^0;

ub1 = u(cv1);
ub2 = u(cv2);
ub3 = u(cv3);
ub4 = u(cv4);
ubi = u(i);
uL = obj.bcRightVal;
if(p>2)
 b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
b = (b-[(A(1,1)/(xbi-(h(i)/2)))*(u(i)-uL) ; 
       (A(2,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
       (A(3,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
       (A(4,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;]);


   C = zeros(4,p-2);
for jj = 1:p-2
   C(1,jj) = (A(1,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1)); 
   C(2,jj) = (A(2,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
   C(3,jj) = (A(3,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
   C(4,jj) = (A(4,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
end
A = (A(:,2:p-1)-C);
            

y(3:p) = (A'*A)\(A'*b);
end
P = [ 1 h(i)/2; 1 xbi];
q = [uL; ubi];

for j = 2:p-1
q = q-[ (h(i)/2)^j*y(j+1) ; obj.moments(i,j+1)*y(j+1)];
end





y(1:2) = P\q;

Z(:,i) = y;





    for i = 3:N

switch i
      case 3
               cv1 =i-1;
cv2 = i+1;
cv3 = i+2;
cv4 = i+3;

    case N
                cv1 =i-1;
cv2 = i+1;
cv3 = i-2;
cv4 = i-3;
  
    

    otherwise
                cv1 =i-1;
cv2 = i+1;
cv3 = i-2;
cv4 = i+2;
end
wi1 = 1/abs(x(cv1)-x(i))^0;
wi2 = 1/abs(x(cv2)-x(i))^0;
wi3 = 1/abs(x(cv3)-x(i))^0;
wi4 = 1/abs(x(cv4)-x(i))^0;

ub1 = u(cv1);
ub2 = u(cv2);
ub3 = u(cv3);
ub4 = u(cv4);
ubi = u(i);


b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];

% AD = obj.primalPI;

Y(2:p) = AD(:,:,i)*b;




 Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
%q = y(1)-ubi

for k = 1:p-1
   Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
end


 Z(:,i) = Y;
 
 
 
    end 
end
end



