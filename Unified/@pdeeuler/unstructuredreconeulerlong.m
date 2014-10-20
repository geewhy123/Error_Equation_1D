function [ Z ] = unstructuredreconeulerlong(obj,u,p,eqn,iUnk)
%UNSTRUCTUREDRECONEULER Summary of this function goes here
%   Detailed explanation goes here

%   [ Z ]= higherunstructuredreconeuler(obj,u,p,eqn,iUnk) ;
% return;

    if(strcmp(eqn,'error')==1)
        [ Z ]= unstructuredreconerroreulerlong(obj,u,p,eqn,iUnk) ;
        return;
    end

%UNSTRUCTUREDRECON3 Summary of this function goes here
%   Detailed explanation goes here
    N = obj.nCells;
    x = obj.cellCentroids;
    h = obj.cellWidths;
    moments = obj.moments;
    Z = zeros(p,N+2);
    wt = obj.weight;


% elseif(obj.bcLeftType=='D' && obj.bcRightType == 'D')
%     
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
    % uL = obj.bcLeftVal(iUnk);


    if(iUnk ~=  3 )
        gam = 1.4;
        Pa = obj.bcLeftVal(3);


% error('1')
        T0=obj.T0;
        P0=obj.P0;
        Ta = T0*(Pa/P0)^((gam-1)/gam);
        rhoa = Pa/Ta;%*gam
        ua = sqrt((2/(gam-1))*(T0/Ta-1)) *sqrt(gam*Pa/rhoa);%???
% % %     fprintf('check quant definitions, and consistent')

        if(strcmp(eqn,'error')==1 && T0 ==0 && P0 == 0)
% [P0 T0 rhoa ua Ta Pa]
            rhoa = 0;
            ua = 0;
            Ta = 0;
            Pa = 0;
            % [P0 T0 rhoa ua Ta Pa]
            %      error('1')
        end

        if(iUnk==1)
            obj.bcLeftVal(1) = rhoa;
            uL = rhoa;
        elseif(iUnk==2)
            obj.bcLeftVal(2) = ua;
            uL = ua;
%     else
%         uL = 0.95;
        end

        A = AA(:,:,i);
        xbi = obj.moments(i,2);

        if(p>2)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi); wi5*(ub5-ubi); wi6*(ub6-ubi) ];

            b = (b-[(A(1,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ; 
            (A(2,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(3,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(4,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(5,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;
            (A(6,1)/(xbi-(-h(i)/2)))*(u(i)-uL) ;]);

            C = zeros(4,p-2);
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

    else 
        b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi5*(ub5-ubi);wi6*(ub6-ubi) ];

% AD = obj.primalPI;
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
   
        PRB = 0;
        for k = 1:p
            PRB = PRB+ Z(k,i)*(-h(i)/2)^(k-1); 
        end
        obj.bcLeftVal(3) = PRB;
 
%  PRB
%  Z
%  obj.bcLeftVal
%  error('1')
    end

% % % 

% %tmp
% b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% 
% % AD = obj.primalPI;
% 
% Y(2:p) = AD(:,:,i)*b;
%  Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
% %q = y(1)-ubi
% 
% for k = 1:p-1
%    Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
% end
% 
% 
%  Z(:,i) = Y;
%%%


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
% % % uL = obj.bcRightVal(iUnk);

% % % 
    if(iUnk == 3)
%     obj.bcRightVal(3) = obj.Pb;
        uL = obj.Pb;%0.97;obj.bcRightVal(3);

        if(p>2)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi);wi4*(ub5-ubi);wi4*(ub6-ubi) ];
            b = (b-[(A(1,1)/(xbi-(h(i)/2)))*(u(i)-uL) ; 
            (A(2,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(3,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(4,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(5,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;
            (A(6,1)/(xbi-(h(i)/2)))*(u(i)-uL) ;]);

            C = zeros(4,p-2);
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

    else
        b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ;wi5*(ub5-ubi);wi6*(ub6-ubi)];
 
        if(obj.stencilSize == 3)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); ];
        end
            
        Y(2:p) = AD(:,:,i)*b;

        Y(1) = ubi;

        for k = 1:p-1
            Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
        end

        Z(:,i) = Y;
  
        if(iUnk == 1)
            rhoRB = 0;
            for k = 1:p
                rhoRB = rhoRB+ Z(k,i)*(h(i)/2)^(k-1); 
            end
            obj.bcRightVal(1) = rhoRB;
   
        elseif(iUnk==2)
     
            uRB = 0;
            for k = 1:p
                uRB = uRB+ Z(k,i)*(h(i)/2)^(k-1); 
            end
            obj.bcRightVal(2) = uRB;
 
        end
        
    end
%%%

% %tmp
% b = [wi1*(ub1-ubi); wi2*(ub2-ubi); wi3*(ub3-ubi); wi4*(ub4-ubi) ];
% 
% % AD = obj.primalPI;
% 
% Y(2:p) = AD(:,:,i)*b;
% 
% 
% 
% 
%  Y(1) = ubi;%-xbi*y(2);%ubi-xbi*y(2)
% %q = y(1)-ubi
% 
% for k = 1:p-1
%    Y(1) = Y(1) - Y(k+1)*moments(i,k+1); 
% end
% 
% 
%  Z(:,i) = Y;
%%%

    
    
    if(strcmp(obj.bchandle,'HC')~=1)
        start = 2;
        fin = N+1;
    else
        start = 3;
       fin = N;
    end

    for i = start:fin

        switch i
            case 2
                cv1 = i+1;
                cv2 = i+2;
                cv3 = i+3;
                cv4 = i+4;
                cv5 = i+5;
                cv6 = i+6;
            case N+1
                cv1 = i-1;
                cv2 = i-2;
                cv3 = i-3; 
                cv4 = i-4;
                cv5 = i-5;
                cv6 = i-6;
              
        case 3
            cv1 =i-1;
            cv2 = i+1;
            cv3 = i+2;
            cv4 = i+3;
            cv5 = i+4;
            cv6 = i+5;

        case N
            cv1 =i-1;
            cv2 = i+1;
            cv3 = i-2;
            cv4 = i-3;
            cv5 = i-4;
            cv6 = i-5;
            case 4
            cv1 =i-1;
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

        if(obj.stencilSize == 3)
            b = [wi1*(ub1-ubi); wi2*(ub2-ubi); ];
        end

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

