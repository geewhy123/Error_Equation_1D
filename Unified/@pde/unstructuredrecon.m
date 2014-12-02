function [ Z] = unstructuredrecon(obj,u,p,eqn)
%UNSTRUCTUREDRECON3 Summary of this function goes here
    if(strcmp(obj.physics,'LinearSystem')==1)
            [Z1] = unstructuredreconlinearsystem (obj,u(:,1),p,eqn,1); 
            [Z2] = unstructuredreconlinearsystem (obj,u(:,2),p,eqn,2);
            Z = [Z1; Z2];
            return;
 
        
    elseif(strcmp(obj.physics,'EulerQ')==1)
        if( (strcmp(eqn,'solution')==1 && obj.hOrder > obj.pOrder) || (strcmp(eqn,'residual')==1 && obj.hOrder > obj.rOrder) || (strcmp(eqn,'error')==1 && obj.hOrder > obj.qOrder) )
            [Z3] = higherunstructuredreconeuler (obj,u(:,3),obj.hOrder,eqn,3);                          
            [Z1] = higherunstructuredreconeuler (obj,u(:,1),obj.hOrder,eqn,1); 
            [Z2] = higherunstructuredreconeuler (obj,u(:,2),obj.hOrder,eqn,2);
            Z = [Z1; Z2;Z3];
            return;
        end
        
        if(p<6)
            [Z3] = unstructuredreconeuler (obj,u(:,3),p,eqn,3);                          
            [Z1] = unstructuredreconeuler (obj,u(:,1),p,eqn,1); 
            [Z2] = unstructuredreconeuler (obj,u(:,2),p,eqn,2);
            Z = [Z1; Z2;Z3];
            return;
        elseif(p==6)
            [Z3] = unstructuredreconeulerlong (obj,u(:,3),p,eqn,3);                          
            [Z1] = unstructuredreconeulerlong (obj,u(:,1),p,eqn,1); 
            [Z2] = unstructuredreconeulerlong (obj,u(:,2),p,eqn,2);
            Z = [Z1; Z2;Z3];
            return;
        else
%                 [Z3] = unstructuredreconeulerlong (obj,u(:,3),p,eqn,3);                          
%                [Z1] = unstructuredreconeulerlong (obj,u(:,1),p,eqn,1); 
%                 [Z2] = unstructuredreconeulerlong (obj,u(:,2),p,eqn,2);
%                Z = [Z1; Z2;Z3];

        return ;
        end
            
    end
    
%  if( (strcmp(eqn,'solution')==1 && obj.hOrder > obj.pOrder) || (strcmp(eqn,'residual')==1 && obj.hOrder > obj.rOrder) || (strcmp(eqn,'error')==1 && obj.hOrder > obj.qOrder) )
%       [Z] = higherunstructuredreconp (obj,u,obj.hOrder,eqn); 
% return;
%  end

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
%    [~,Z] = unstructuredrecon5 (u,x,h,N,u0,u1); 
        [Z] = unstructuredreconlong (obj,u,p,eqn); 
    otherwise 
        
        [Z] = unstructuredreconp (obj,u,p,eqn); 
    
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
    wt = obj.weight;

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

            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;
    
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
            if(obj.hOrder > 0)
                AA = obj.higherprimalRM;
                AD = obj.higherprimalPI;
            end
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
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;

        ub1 = u(cv1);
        ub2 = u(cv2);
        ub3 = u(cv3);
        ub4 = u(cv4);
        ubi = u(i);
        uL = obj.bcLeftVal;

        A = AA(:,:,i);
        xbi = obj.moments(i,2);

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


% [m,n] = size(A);
% if(n == 1)
%     i
%     AA
%     eqn
%     obj.errorRM
% end

            A = (A(:,2:p-1)-C);
        

            y(3:p) = (A'*A)\(A'*b);
        end

        P = [ 1 -h(i)/2; 1 xbi];
        q = [uL; ubi];
        for j = 2:p-1
            q = q-[ (-h(i)/2)^j*y(j+1) ; obj.moments(i,j+1)*y(j+1)];
        end

        y(1:2) = P\q;

%%%%%
% 
% A(2:5) = A(1:4);
% A(1) = -h(i)/2-xbi;
%  b = [wi1*(uL-ubi); wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% 
% y(2) = (A'*A)\(A'*b);
% y(1) = ubi-y(2)*xbi;
% % y
% % error('1')
%%%%%


        Z(:,i) = y;


        i = N+1;
        xbi = obj.moments(i,2);

        A = AA(:,:,i);
        cv1 =i-1;
        cv2 = i-2;
        cv3 = i-3;
        cv4 = i-4;
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;

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


%%%%%
% 
% A(2:5) = A(1:4);
% A(1) = h(i)/2-xbi;
%  b = [wi1*(uL-ubi); wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% 
% y(2) = (A'*A)\(A'*b);
% y(1) = ubi-y(2)*xbi;
% % y
% % error('1')
% %%%%%
        Z(:,i) = y;

        if(strcmp(obj.bchandle,'HC')==1)
           begin = 3;
           fin = N;
        elseif(strcmp(obj.bchandle,'Jump')==1)
           begin = 2;
           fin = N+1;
        end
               
        for i = begin:fin    
            switch i
            case 2
                cv1 = i+1;
                cv2 = i+2;
                cv3 = i+3;
                cv4 = i+4;
            case N+1
                cv1 = i-1;
                cv2 = i-2;
                cv3 = i-3;
                cv4 = i-4;
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
            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;

            ub1 = u(cv1);
            ub2 = u(cv2);
            ub3 = u(cv3);
            ub4 = u(cv4);
            ubi = u(i);

            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];

% AD = obj.primalPI;
% size(AD(:,:,i)*b)
% p
% i
            if(obj.stencilSize == 3)
                b = [wi1*(ub1-ubi); wi2*(ub2-ubi); ];
            end

            Y(2:p) = AD(:,:,i)*b;

            Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
%q = y(1)-ubi

            for k = 1:p-1
                Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
            end

            Z(:,i) = Y;
 
        end
    
    elseif(obj.bcLeftType=='F' && obj.bcRightType == 'D')
    
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
    
    
% % % i=2;
% % % cv1 =i+1;
% % % cv2 = i+2;
% % % cv3 = i+3;
% % % cv4 = i+4;
% % % wi1 = 1/abs(x(cv1)-x(i))^wt;
% % % wi2 = 1/abs(x(cv2)-x(i))^wt;
% % % wi3 = 1/abs(x(cv3)-x(i))^wt;
% % % wi4 = 1/abs(x(cv4)-x(i))^wt;
% % % 
% % % ub1 = u(cv1);
% % % ub2 = u(cv2);
% % % ub3 = u(cv3);
% % % ub4 = u(cv4);
% % % ubi = u(i);
% % % uL = obj.bcLeftVal;
% % % 
% % % 
% % % 
% % % A = AA(:,:,i);
% % % xbi = obj.moments(i,2);
% % % x2bi = obj.moments(i,3);
% % % if(p>2)
% % %  b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% % %  
% % % b = (b-[(A(1,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ; 
% % %        (A(2,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
% % %        (A(3,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
% % %        (A(4,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;]);
% % %     
% % % %    
% % % % C =           ([(A(1,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2) ; 
% % % %                 (A(2,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2)  ;
% % % %                 (A(3,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2)  ;
% % % %                 (A(4,1)/(xbi-(-h(i)/2)))*(x2bi-(-h(i)/2)^2) ]);   
% % % C = zeros(4,p-2);
% % % for jj = 1:p-2
% % %    C(1,jj) = (A(1,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1)); 
% % %    C(2,jj) = (A(2,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
% % %    C(3,jj) = (A(3,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
% % %    C(4,jj) = (A(4,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
% % % end
% % % 
% % % A = (A(:,2:p-1)-C);
% % %         
% % % 
% % % y(3:p) = (A'*A)\(A'*b);
% % % end
% % % 
% % % P = [ 1 -h(i)/2; 1 xbi];
% % % q = [uL; ubi];
% % % for j = 2:p-1
% % % q = q-[ (-h(i)/2)^j*y(j+1) ; obj.moments(i,j+1)*y(j+1)];
% % % end
% % % 
% % % 
% % % 
% % % 
% % % y(1:2) = P\q;
% % % 
% % % Z(:,i) = y;


        i = N+1;
        xbi = obj.moments(i,2);
        A = AA(:,:,i);
        cv1 =i-1;
        cv2 = i-2;
        cv3 = i-3;
        cv4 = i-4;
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;

        ub1 = u(cv1);
        ub2 = u(cv2);
        ub3 = u(cv3);
        ub4 = u(cv4);
        ubi = u(i);

        obj.curTime;
        uL = obj.bcRightVal;%sin(pi*obj.curTime);%obj.bcRightVal;%-1/exp(obj.step*obj.tStep)
        obj.bcRightVal = uL;

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
        
        for i = 2:N

            switch i
            case 2
                cv1 =i+1;
                cv2 = i+2;
                cv3 = i+3;
                cv4 = i+4;
  
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
            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;

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



function [Z] = unstructuredreconlong(obj,u ,p,eqn)
%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
    N = obj.nCells;
    x = obj.cellCentroids;
    h = obj.cellWidths;
    moments = obj.moments;
    Z = zeros(p,N+2);
    wt = obj.weight;

    if(obj.bcLeftType == 'P' && obj.bcRightType == 'P')
        if(strcmp(eqn,'solution')==1)
            AD = obj.primalPI;
        elseif(strcmp(eqn,'error')==1)
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
                cv5 = N-1;
                cv6 = i+3;
            case 3
                cv1 = i-1;
                cv2 = i+1;
                cv3 = N+1;
                cv4 = i+2;
                cv5 = N;
                cv6 = i+3;
            case 4
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = N+1;
                cv6 = i+3;
            case N-1
                cv1 = i-1;
                cv2 = N;
                cv3 = i-2;
                cv4 = N+1;
                cv5 = i-3;
                cv6 = 2;
            case N
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = 2;
                cv5 = i-3;
                cv6 = 3;
            case N+1
                cv1 = i-1;
                cv2 = 2;
                cv3 = i-2;
                cv4 = 3;
                cv5 = i-3;
                cv6 = 4;
            otherwise
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i-3;
                cv6 = i+3;
            end

            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;
            wi5 = 1/abs(x(cv5)-x(i))^wt;
            wi6 = 1/abs(x(cv6)-x(i))^wt;
            ub1 = u(cv1);
            ub2 = u(cv2);
            ub3 = u(cv3);
            ub4 = u(cv4);
            ub5 = u(cv5);
            ub6 = u(cv6);
            ubi = u(i);

            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];

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
        if(strcmp(eqn,'solution')==1)
            AA = obj.primalRM;
            AD = obj.primalPI;
        elseif(strcmp(eqn,'error')==1)
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
        cv1 = i+1;
        cv2 = i+2;
        cv3 = i+3;
        cv4 = i+4;
        cv5 = i+5;
        cv6 = i+6;
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;
        wi5 = 1/abs(x(cv5)-x(i))^wt;
        wi6 = 1/abs(x(cv6)-x(i))^wt;
        ub1 = u(cv1);
        ub2 = u(cv2);
        ub3 = u(cv3);
        ub4 = u(cv4);
        ub5 = u(cv5);
        ub6 = u(cv6);
        ubi = u(i);
        uL = obj.bcLeftVal;

        A = AA(:,:,i);
        xbi = obj.moments(i,2);

        if(p>2)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];

            b = (b-[(A(1,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ; 
            (A(2,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(3,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(4,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(5,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(6,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;]);

            C = zeros(6,p-2);
            for jj = 1:p-2
                C(1,jj) = (A(1,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1)); 
                C(2,jj) = (A(2,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
                C(3,jj) = (A(3,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
                C(4,jj) = (A(4,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
                C(5,jj) = (A(5,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
                C(6,jj) = (A(6,1)/(xbi-(-h(i)/2)))*(obj.moments(i,jj+2)-(-h(i)/2)^(jj+1));
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
        A = AA(:,:,i);
        cv1 =i-1;
        cv2 = i-2;
        cv3 = i-3;
        cv4 = i-4;
        cv5 = i-5;
        cv6 = i-6;
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;
        wi5 = 1/abs(x(cv5)-x(i))^wt;
        wi6 = 1/abs(x(cv6)-x(i))^wt;
        ub1 = u(cv1);
        ub2 = u(cv2);
        ub3 = u(cv3);
        ub4 = u(cv4);
        ub5 = u(cv5);
        ub6 = u(cv6);
        ubi = u(i);
        uL = obj.bcRightVal;

        if(p>2)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];
            b = (b-[(A(1,1)/(xbi-(h(i)/2)))*(u(i)-uL) ; 
            (A(2,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(3,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(4,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(5,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(6,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;]);

            C = zeros(6,p-2);

            for jj = 1:p-2
                C(1,jj) = (A(1,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1)); 
                C(2,jj) = (A(2,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(3,jj) = (A(3,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(4,jj) = (A(4,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(5,jj) = (A(5,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(6,jj) = (A(6,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
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
                cv5 = i+4;
                cv6 = i+5;
            case N
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i-3;
                cv5 = i-4;
                cv6 = i-5;
            case 4
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i+3;
                cv6 = i+4;        
            case N-1
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i-3;
                cv6 = i-4;

            otherwise
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i-3;
                cv6 = i+3;
            end
            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;
            wi5 = 1/abs(x(cv5)-x(i))^wt;
            wi6 = 1/abs(x(cv6)-x(i))^wt;
            ub1 = u(cv1);
            ub2 = u(cv2);
            ub3 = u(cv3);
            ub4 = u(cv4);
            ub5 = u(cv5);
            ub6 = u(cv6);
            ubi = u(i);


            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];

% AD = obj.primalPI;

            Y(2:p) = AD(:,:,i)*b;
            Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
%q = y(1)-ubi

            for k = 1:p-1
                Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
            end


            Z(:,i) = Y;

        end

    elseif(obj.bcLeftType=='F' && obj.bcRightType == 'D')

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


        i = N+1;
        xbi = obj.moments(i,2);
        A = AA(:,:,i);
        cv1 =i-1;
        cv2 = i-2;
        cv3 = i-3;
        cv4 = i-4;
        cv5 = i-5;
        cv6 = i-6;
        wi1 = 1/abs(x(cv1)-x(i))^wt;
        wi2 = 1/abs(x(cv2)-x(i))^wt;
        wi3 = 1/abs(x(cv3)-x(i))^wt;
        wi4 = 1/abs(x(cv4)-x(i))^wt;
        wi5 = 1/abs(x(cv5)-x(i))^wt;
        wi6 = 1/abs(x(cv6)-x(i))^wt;

        ub1 = u(cv1);
        ub2 = u(cv2);
        ub3 = u(cv3);
        ub4 = u(cv4);
        ub5 = u(cv5);
        ub6 = u(cv6);
        ubi = u(i);

        uL = obj.bcRightVal;%sin(pi*obj.curTime);%obj.bcRightVal;%-1/exp(obj.step*obj.tStep)
        obj.bcRightVal = uL;

        if(p>2)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];
            b = (b-[(A(1,1)/(xbi-(h(i)/2)))*(u(i)-uL) ; 
            (A(2,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(3,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(4,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(5,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(6,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;]);

            C = zeros(6,p-2);
            for jj = 1:p-2
                C(1,jj) = (A(1,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1)); 
                C(2,jj) = (A(2,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(3,jj) = (A(3,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(4,jj) = (A(4,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(5,jj) = (A(5,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
                C(6,jj) = (A(6,1)/(xbi-(h(i)/2)))*(obj.moments(i,jj+2)-(h(i)/2)^(jj+1));
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


        for i = 2:N

            switch i
            case 2
                cv1 =i+1;
                cv2 = i+2;
                cv3 = i+3;
                cv4 = i+4;
                cv5 = i+5;
                cv6 = i+6;


            case 3
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i+2;
                cv4 = i+3;
                cv5 = i+4;
                cv6 = i+5;
            case 4
                cv1 = i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i+3;
                cv6 = i+4;

            case N-1
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i-3;
                cv6 = i-4;        
            case N
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i-3;
                cv5 = i-4;
                cv6 = i-5;

            otherwise
                cv1 =i-1;
                cv2 = i+1;
                cv3 = i-2;
                cv4 = i+2;
                cv5 = i-3;
                cv6 = i+3;
            end
            wi1 = 1/abs(x(cv1)-x(i))^wt;
            wi2 = 1/abs(x(cv2)-x(i))^wt;
            wi3 = 1/abs(x(cv3)-x(i))^wt;
            wi4 = 1/abs(x(cv4)-x(i))^wt;
            wi5 = 1/abs(x(cv5)-x(i))^wt;
            wi6 = 1/abs(x(cv6)-x(i))^wt;
            ub1 = u(cv1);
            ub2 = u(cv2);
            ub3 = u(cv3);
            ub4 = u(cv4);
            ub5 = u(cv5);
            ub6 = u(cv6);
            ubi = u(i);

            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi); wi5*(ub5-ubi); wi6*(ub6-ubi) ];

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



