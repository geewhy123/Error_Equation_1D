classdef pde < handle
    %PDE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nCells;
        pOrder;
        qOrder;
         rOrder;
         bcLeftType;
         bcLeftVal;
         bcRightType;
         bcRightVal;
         endTime;
         tOrder;
         physics;
         goal;
         cellCentroids;
         cellWidths;
         initialSolution;
         exactSolution;
         source;
         primalPI;
         moments; 
         primalRM;
         xhat;
         tStep;
         resPI;
         resRM;
         errorPI;
         errorRM;
         residual;
         Rsp;
         errorSource;
         curTime;
         weight;
         convSoln;
         convSolnRecon;
higherprimalPI;
        higherprimalRM;
        hOrder=0;
        refinecells;
    end
    
    methods
        function obj = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k,wt)
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
%     end
%     methods( Static=true)
        initializeexact(obj);
        computemoments(obj);
        computeprimalpseudo(obj);
        Z = unstructuredrecon(obj,u,order,eqn);
        er = reconplot(obj,Z,eqn);
        [uu,d] = updatesolution(obj,u);
        FI = computefluxintegral(obj,u,eqn);
        J=computefluxjacobian(obj,u,eqn);
        R = computeres(obj,u,time,r);
        computerespseudo(obj);
        [AD,AA] = computepseudo(obj,p);
        computeerrorpseudo(obj);
        [errerr2,x,cverr2,exacterr,ee,te  ]=solvebyjacobian(obj);
        [errerr2,x,cverr2,exacterr,ee,te  ]=solvebyjacobianNL(obj);
        computehigherpseudo(obj);
    end
    
    methods(Access=private)
       
        poissoninitialize(obj);
        advectioninitialize(obj);
        burgersinitialize(obj);
    end
        
    
end

