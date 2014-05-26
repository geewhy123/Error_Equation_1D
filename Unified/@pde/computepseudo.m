function  [AD,AA] = computepseudo(obj,p)
%COMPUTEPSEUDO Summary of this function goes here
%   Detailed explanation goes here
N = obj.nCells;
h = obj.cellWidths;
x = obj.cellCentroids;

% p = obj.pOrder;
AD = zeros(p-1,4,N+2);
AA = zeros(4,p-1,N+2);
for i = 2:N+1
    x1 = 0;
    x2 = 0;
    x3 = 0;
    x4 = 0;
   
    
    if ((i > 3) && (i < N))
              cv1 = i-1;
        cv2 = i+1;
        cv3 = i-2;
        cv4 = i+2;
  

    end

  if i==3
            cv1 = i-1;
        cv2 = i+1;
        cv3 = N+1;
        x3 = x3-1;
        cv4 = i+2;

    if(obj.bcRightType == 'D')
        x3 = 0;
        cv3 = i+2;
        cv4 = i+3;
   
    end
   end
    if i==N
        cv1 = i-1;
        cv2 = i+1;
        cv3 = i-2;
        cv4 = 2;
        x4 = x4+1;

     if(obj.bcRightType == 'D')
         x4 = 0;
         cv4 = i-3;
    end
    end
    if i==2
        cv1 = N+1;
        x1 = x1 -1;
        cv2 = i+1;
        cv3 = N;
        x3 = x3-1;
        cv4 = i+2;

    if(obj.bcLeftType == 'D' || obj.bcLeftType == 'F')
        x1 = 0;
        x3 = 0;
        cv1 = i+1;
        cv2 = i+2;
        cv3 = i+3;
        cv4 = i+4;
    else
%          assert(0)
        
    end
    
    
    
    end
    if i==N+1
        cv1 = i-1;
        cv2 = 2;
        x2 = x2 +1;
        cv3 = i-2;
        cv4 = 3;
        x4 = x4+1;

        if(obj.bcRightType == 'D')
            x2 = 0;
            x4 = 0;
            
            cv1 = i-1;
        cv2 = i-2;
        cv3 = i-3;
        cv4 = i-4;
    end
    end
    
    x1 = x1+x(cv1);
    x2 = x2+x(cv2);
    x3 = x3+x(cv3);
    x4 = x4+x(cv4);
    xi = x(i);
   
wi1 = 1/abs(x1-xi)^0;
wi2 = 1/abs(x2-xi)^0;
wi3 = 1/abs(x3-xi)^0;
wi4 = 1/abs(x4-xi)^0;
    

% x2b1

 
 xhat = zeros(4,p-1);
for k = 1:p-1
   for ii = 1:k+1
%       if(k==2 && ii==1)
%           i
%          [x2b1 nchoosek(k,ii-1)*(x1-xi)^(ii-1)*obj.moments(cv1,k-ii+2)]
%       end
xhat(1,k) = xhat(1,k)+ nchoosek(k,ii-1)*(x1-xi)^(ii-1)*obj.moments(cv1,k-ii+2);   
xhat(2,k) = xhat(2,k)+nchoosek(k,ii-1)*(x2-xi)^(ii-1)*obj.moments(cv2,k-ii+2);
xhat(3,k) = xhat(3,k)+nchoosek(k,ii-1)*(x3-xi)^(ii-1)*obj.moments(cv3,k-ii+2);
xhat(4,k) = xhat(4,k)+nchoosek(k,ii-1)*(x4-xi)^(ii-1)*obj.moments(cv4,k-ii+2);
   end

        
end
%  [x2b1+2*(x1-xi)*xb1+(x1-xi)^2 c]

   
   

for j = 1:p-1
   AA(1,j,i) = xhat(1,j)-obj.moments(i,j+1); 
   AA(2,j,i) = xhat(2,j)-obj.moments(i,j+1);
   AA(3,j,i) = xhat(3,j)-obj.moments(i,j+1);
   AA(4,j,i) = xhat(4,j)-obj.moments(i,j+1);
end

% AA;
% BB
% AA = double([wi1*(xb1-xbi+x1-xi) wi1*(x2b1+2*(x1-xi)*xb1+(x1-xi)^2-x2bi) ; 
%              wi2*(xb2-xbi+x2-xi) wi2*(x2b2+2*(x2-xi)*xb2+(x2-xi)^2-x2bi) ; 
%              wi3*(xb3-xbi+x3-xi) wi3*(x2b3+2*(x3-xi)*xb3+(x3-xi)^2-x2bi) ;
%              wi4*(xb4-xbi+x4-xi) wi4*(x2b4+2*(x4-xi)*xb4+(x4-xi)^2-x2bi)  ])

 AD(:,:,i)=pinv(AA(:,:,i));
end


% AD
% assert(0)

%  error('1')
% obj.primalPI = AD;
% obj.reconM = AA;
obj.xhat = xhat;
end

