function [ phi1,phi2,phi3 ] = computeeulerdifffluxintegral4Newton( obj,Z,eqn )
%COMPUTEEULERFLUXINTEGRAL Summary of this function goes here
%   Detailed explanation goes here
% error('2')
    N=obj.nCells;
    h = obj.cellWidths;
    x = obj.cellCentroids;
    rhol  = zeros(N+2,1);
    rhor  = zeros(N+2,1);
    ul = zeros(N+2,1);
    ur = zeros(N+2,1);


    Pl = zeros(N+2,1);
    Pr = zeros(N+2,1);
    rhob= zeros(N+2,1);
    ub= zeros(N+2,1);
    Pb= zeros(N+2,1);

    Ut1 = zeros(N+2,1);
    Ut2 = zeros(N+2,1);
    Ut3 =zeros(N+2,1);
    gam = 1.4;
    bAtilde = zeros(3,3,N+2);
    FrAve = zeros(N+2,3);
    FlAve = zeros(N+2,3);

    if(strcmp(eqn,'solution')==1)
        order = obj.pOrder; 
    elseif(strcmp(eqn,'residual')==1)
        order = obj.rOrder;
    elseif(strcmp(eqn,'error')==1)
        order = obj.qOrder;
    else
        assert(0)
    end


    for i = 2:N+1
    
        for k = 1:order
            rhor(i) = rhor(i)+ Z(k,i)*(h(i)/2)^(k-1);
            rhol(i) = rhol(i)+ Z(k,i)*(-h(i)/2)^(k-1);
            ur(i)   = ur(i)+ Z(k+order,i)*(h(i)/2)^(k-1);
            ul(i)   = ul(i)+ Z(k+order,i)*(-h(i)/2)^(k-1);
            Pr(i)   = Pr(i)+ Z(k+2*order,i)*(h(i)/2)^(k-1);
            Pl(i)   = Pl(i)+ Z(k+2*order,i)*(-h(i)/2)^(k-1);

            rhob(i) = rhob(i) + (Z(k,i))*obj.moments(i,k);
            ub(i) = ub(i) + (Z(k+order,i))*obj.moments(i,k);
            Pb(i) = Pb(i) + (Z(k+2*order,i))*obj.moments(i,k);
        end


    end




    U1l = obj.convUleft(:,1);
    U2l = obj.convUleft(:,2);
    U3l = obj.convUleft(:,3);
    U1r = obj.convUright(:,1);
    U2r = obj.convUright(:,2);
    U3r = obj.convUright(:,3);
        
    F1l = zeros(N+2,1);
    F2l = zeros(N+2,1);
    F3l = zeros(N+2,1);
    F1r = zeros(N+2,1);
    F2r = zeros(N+2,1);
    F3r = zeros(N+2,1);
Ql  = zeros(3,3);
Qr = zeros(3,3);
    for i = 2:N+1
        Ql(1,1) = 0;
        Ql(1,2) = 1;
        Ql(1,3) = 0;
        Ql(2,1) = ((gam-1)/2-1)*U2l(i)^2/U1l(i)^2;
        Ql(2,2) = (3-gam)*U2l(i)/U1l(i);
        Ql(2,3) = (gam-1);
        Ql(3,1) = (gam-1)*U2l(i)^3/U1l(i)^3-gam*U2l(i)*U3l(i)/U1l(i)^2;
        Ql(3,2) = gam*U3l(i)/U1l(i) - 1.5*(gam-1)*U2l(i)^2/U1l(i)^2;
        Ql(3,3) = gam*U2l(i)/U1l(i);
        
        Qr(1,1) = 0;
        Qr(1,2) = 1;
        Qr(1,3) = 0;
        Qr(2,1) = ((gam-1)/2-1)*U2r(i)^2/U1r(i)^2;
        Qr(2,2) = (3-gam)*U2r(i)/U1r(i);
        Qr(2,3) = (gam-1);
        Qr(3,1) = (gam-1)*U2r(i)^3/U1r(i)^3-gam*U2r(i)*U3r(i)/U1r(i)^2;
        Qr(3,2) = gam*U3r(i)/U1r(i) - 1.5*(gam-1)*U2r(i)^2/U1r(i)^2;
        Qr(3,3) = gam*U2r(i)/U1r(i);
        
        Fl = Ql*[rhol(i) ul(i) Pl(i)]';
        F1l(i) = Fl(1);
        F2l(i) = Fl(2);
        F3l(i) = Fl(3);
        
        Fr = Qr*[rhor(i) ur(i) Pr(i)]';
        F1r(i) = Fr(1);
        F2r(i) = Fr(2);
        F3r(i) = Fr(3);
    end
    
    
%         [F1l(2),F2l(2),F3l(2)]=conservedtoflux(U1l(2),U2l(2),U3l(2));
%         [F1r(2),F2r(2),F3r(2)]=conservedtoflux(U1r(2),U2r(2),U3r(2));
%         [F1l(N+1),F2l(N+1),F3l(N+1)]=conservedtoflux(U1l(N+1),U2l(N+1),U3l(N+1));
%         [F1r(N+1),F2r(N+1),F3r(N+1)]=conservedtoflux(U1r(N+1),U2r(N+1),U3r(N+1));
%
%
%     for i = 2:N+1
%         [F1l(i),F2l(i),F3l(i)]=conservedtoflux(U1l(i),U2l(i),U3l(i));
%         [F1r(i),F2r(i),F3r(i)]=conservedtoflux(U1r(i),U2r(i),U3r(i));
%     end

    if(strcmp(obj.bchandle,'HC')~=1)
        obj.P0;
        obj.T0;
        obj.Pb;
% [F1l(2),F2l(2),F3l(2)]=inboundaryflux(U1l(2)+rhol(2),U2l(2)+ul(2),U3l(2)+Pl(2),obj.primalP0,obj.primalT0);
% [F1r(N+1),F2r(N+1),F3r(N+1)]=outboundaryflux(U1r(N+1)+rhor(N+1),U2r(N+1)+ur(N+1),U3r(N+1)+Pr(N+1),obj.primalPb);
% [F1lu(2),F2lu(2),F3lu(2)]=inboundaryflux(U1l(2),U2l(2),U3l(2),obj.primalP0,obj.primalT0);
% [F1ru(N+1),F2ru(N+1),F3ru(N+1)]=outboundaryflux(U1r(N+1),U2r(N+1),U3r(N+1),obj.primalPb);
% F1l(2) = F1l(2) - F1lu(2);
% F2l(2) = F2l(2) - F2lu(2);
% F3l(2) = F3l(2) - F3lu(2);
% F1r(N+1) = F1r(N+1) - F1ru(N+1);
% F2r(N+1) = F2r(N+1) - F2ru(N+1);
% F3r(N+1) = F3r(N+1) - F3ru(N+1);
[F1l(2),F2l(2),F3l(2)]=inboundaryfluxNewton(U1l(2),U2l(2),U3l(2),rhol(2),ul(2),Pl(2),obj.primalP0,obj.primalT0);
[F1r(N+1),F2r(N+1),F3r(N+1)]=outboundaryfluxNewton(U1r(N+1),U2r(N+1),U3r(N+1),rhor(N+1),ur(N+1),Pr(N+1),obj.primalPb);


end


    for i = 2:N
        [Ut1(i),Ut2(i),Ut3(i)]=computeroeavg(U1l(i+1),U2l(i+1),U3l(i+1),U1r(i),U2r(i),U3r(i));    

        bAtilde(:,:,i) = computeAtilde(Ut1(i),Ut2(i),Ut3(i));
    end


    for i = 2:N+1
        FrAve(i,1:3) =(0.5*[(F1r(i)+F1l(i+1)); (F2r(i)+F2l(i+1)); (F3r(i)+F3l(i+1))]  -0.5*bAtilde(:,:,i)*  ([ rhol(i+1); ul(i+1); Pl(i+1)]- [ rhor(i); ur(i); Pr(i)]) )';
        FlAve(i,1:3) =(0.5*[(F1l(i)+F1r(i-1)); (F2l(i)+F2r(i-1)); (F3l(i)+F3r(i-1))]  -0.5*bAtilde(:,:,i-1)*([ rhol(i); ul(i); Pl(i)]- [ rhor(i-1); ur(i-1); Pr(i-1)]) )';
    end

    FlAve(2,1:3) = [F1l(2); F2l(2); F3l(2)]';
    FrAve(N+1,1:3) = [F1r(N+1);F2r(N+1);F3r(N+1)]';

    F = [FlAve FrAve];

    PAp = NaN*ones(N+2,1);

    A = zeros(N+2,1);
    phib1 = zeros(N+2,1);
    phib2 = zeros(N+2,1);
    phib3 = zeros(N+2,1);

    c1 = 0.3478548451;
    c2 = 0.6521451549;
    c3 = 0.6521451549;
    c4 = 0.3478548451;
    x1= 0.8611363116;
    x2 = 0.339981436;
    x3 = -0.339981436;
    x4= -0.8611363116;


    for i = 2:N+1
        xr = x(i)+h(i)/2;   
        A(i)=obj.getArea(xr);
    end
    A(1) = obj.getArea(0);
 
 ZE = Z;
 Z = obj.unstructuredrecon(obj.convSoln,order,'error');
 
 S = zeros(3,3);
    for i = 2:N+1
        xl = x(i)-h(i)/2;
        xr = x(i)+h(i)/2;
        
        xx1 = ((xr-xl)/2)*x1+(xr+xl)/2;
        xx2 = ((xr-xl)/2)*x2+(xr+xl)/2;
        xx3 = ((xr-xl)/2)*x3+(xr+xl)/2;
        xx4 = ((xr-xl)/2)*x4+(xr+xl)/2;
  
        r1 = 0;
        r2 = 0;
        r3 = 0;
        r4 = 0;
        ru1 = 0;
        ru2 = 0;
        ru3 = 0;
        ru4 = 0;
        rE1 = 0;
        rE2 = 0;
        rE3 = 0;
        rE4 = 0;
        er1 = 0;
        er2 = 0;
        er3 = 0;
        er4 = 0;
        eru1 = 0;
        eru2 = 0;
        eru3 = 0;
        eru4 = 0;
        erE1 = 0;
        erE2 = 0;
        erE3 = 0;
        erE4 = 0;
        for k = 1:order
            r1 = r1 + Z(k,i)*(xx1-x(i))^(k-1); 
            r2 = r2 + Z(k,i)*(xx2-x(i))^(k-1) ;
            r3 = r3 + Z(k,i)*(xx3-x(i))^(k-1) ;
            r4 = r4 + Z(k,i)*(xx4-x(i))^(k-1) ;
            ru1 = ru1 + Z(k+order,i)*(xx1-x(i))^(k-1); 
            ru2 = ru2 + Z(k+order,i)*(xx2-x(i))^(k-1) ;
            ru3 = ru3 + Z(k+order,i)*(xx3-x(i))^(k-1) ;
            ru4 = ru4 + Z(k+order,i)*(xx4-x(i))^(k-1) ;
            rE1 = rE1 + Z(k+2*order,i)*(xx1-x(i))^(k-1); 
            rE2 = rE2 + Z(k+2*order,i)*(xx2-x(i))^(k-1) ;
            rE3 = rE3 + Z(k+2*order,i)*(xx3-x(i))^(k-1) ;
            rE4 = rE4 + Z(k+2*order,i)*(xx4-x(i))^(k-1) ;
            
            er1 = er1 +ZE(k,i)*(xx1-x(i))^(k-1) ;
            er2 = er2 +ZE(k,i)*(xx2-x(i))^(k-1) ;
            er3 = er3 +ZE(k,i)*(xx3-x(i))^(k-1) ;
            er4 = er4 +ZE(k,i)*(xx4-x(i))^(k-1) ;
            eru1 = eru1 +ZE(k+order,i)*(xx1-x(i))^(k-1) ;
            eru2 = eru2 +ZE(k+order,i)*(xx2-x(i))^(k-1) ;
            eru3 = eru3 +ZE(k+order,i)*(xx3-x(i))^(k-1) ;
            eru4 = eru4 +ZE(k+order,i)*(xx4-x(i))^(k-1) ;
            erE1 = erE1 +ZE(k+2*order,i)*(xx1-x(i))^(k-1) ;
            erE2 = erE2 +ZE(k+2*order,i)*(xx2-x(i))^(k-1) ;
            erE3 = erE3 +ZE(k+2*order,i)*(xx3-x(i))^(k-1) ;
            erE4 = erE4 +ZE(k+2*order,i)*(xx4-x(i))^(k-1) ;
            
            
        end
        P1 = (gam-1)*(rE1-0.5*ru1^2/r1);
        P2 = (gam-1)*(rE2-0.5*ru2^2/r2);
        P3 = (gam-1)*(rE3-0.5*ru3^2/r3);
        P4 = (gam-1)*(rE4-0.5*ru4^2/r4);
 
        S(1,:) = 0;
        S(2,1) = ((gam-1)/2)*(ru1^2/r1^2)*obj.getAp(xx1);%/obj.getArea(xx1);
        S(2,2) = -(gam-1)*(ru1/r1)*obj.getAp(xx1);%/obj.getArea(xx1);
        S(2,3) =  (gam-1)*obj.getAp(xx1);%/obj.getArea(xx1);
        S(3,:) = 0;
        
        SE1 = S(2,:)*[er1 eru1 erE1]';
        
        S(2,1) = ((gam-1)/2)*(ru1^2/r1^2)*obj.getAp(xx2);%/obj.getArea(xx2);
        S(2,2) = -(gam-1)*(ru1/r1)*obj.getAp(xx2);%/obj.getArea(xx2);
        S(2,3) =  (gam-1)*obj.getAp(xx2);%/obj.getArea(xx2);
        SE2 = S(2,:)*[er2 eru2 erE2]';
        
        S(2,1) = ((gam-1)/2)*(ru1^2/r1^2)*obj.getAp(xx3);%/obj.getArea(xx3);
        S(2,2) = -(gam-1)*(ru1/r1)*obj.getAp(xx3);%/obj.getArea(xx3);
        S(2,3) =  (gam-1)*obj.getAp(xx3);%/obj.getArea(xx3);
        SE3 = S(2,:)*[er3 eru3 erE3]';
        
        S(2,1) = ((gam-1)/2)*(ru1^2/r1^2)*obj.getAp(xx4);%/obj.getArea(xx4);
        S(2,2) = -(gam-1)*(ru1/r1)*obj.getAp(xx4);%/obj.getArea(xx4);
        S(2,3) =  (gam-1)*obj.getAp(xx4);%/obj.getArea(xx4);
        SE4 = S(2,:)*[er4 eru4 erE4]';
%         [SE1 SE2 SE3 SE4]
        
        PAp(i) = (1/h(i))*(c1*SE1+c2*SE2+c3*SE3+c4*SE4)*(xr-xl)/2;
        phib1(i) =(A(i)*FrAve(i,1)-A(i-1)*FlAve(i,1))/h(i);
        phib2(i) = (A(i)*FrAve(i,2)-A(i-1)*FlAve(i,2))/h(i)- PAp(i);
        phib3(i) =(A(i)*FrAve(i,3)-A(i-1)*FlAve(i,3))/h(i);
    end

% 
%     phi1 = phia1-phib1;
%     phi2 = phia2-phib2;
%     phi3 = phia3-phib3;
    phi1 = phib1;
    phi2 = phib2;
    phi3 = phib3;

    phi1 = phi1 - obj.errorSource(:,1);
    phi2 = phi2 - obj.errorSource(:,2);
    phi3 = phi3 - obj.errorSource(:,3);
% [phi1 phi2 phi3]
%     error('1')

end

