classdef pdeeuler < pde
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gamma = 1.4;
        areatype ;
        P0;
        T0;
        Pb;
        primalP0;
        primalT0;
        primalPb;
        exactSolutionV;
        exactSolutionU;
        convVleft;
        convVright;
        convUleft;
        convUright;
        convVreconp;
        convSolutionV;
        reconexactsolutionV;
%         higherprimalPI;
%         higherprimalRM;
%         hOrder=0;
ebcL;
ebcR;
NLfluxtype=2;
    end
    
    methods
         function obj = pdeeuler(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,wt)
            if(nargin>0)
         obj.nCells = N;
        obj.pOrder = p;
        obj.qOrder = q;
         obj.rOrder = r;
         obj.bcLeftType = BCLeft;
         obj.bcLeftVal = valLeft;
         obj.bcRightType = BCRight;
         obj.bcRightVal = valRight;
         obj.endTime = tlim;
         obj.tOrder = tord;
         obj.physics = physics;
         obj.goal = goal;
         obj.cellCentroids = x;
         obj.cellWidths = h;
         obj.tStep = k;
         obj.weight = wt;
         
            end
         end
        
        
       A = getArea(obj,xx);
       Ap = getAp(obj,xx);
       [xx,rho,u,P] = initializeeuler(obj);
        [errerr2,x,cverr2,exacterr,ee,te  ] =solvebyeulerjacobian(obj);
       [phi1,phi2,phi3]=computeeulerfluxintegral(obj,Z,eqn);
       J = computeeulerfluxjacobian(obj,v,eqn);
        computeprimalleftright(obj);
%         computehigherpseudo(obj);
        

%        [ u1,u2,u3 ] = toconservedvars(obj, rho,u,P );
    end
    
end

