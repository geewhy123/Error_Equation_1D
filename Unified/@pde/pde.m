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
         reconM;
         xhat;
         tStep;
    end
    
    methods
        function obj = pde(N,p,q,r,BCLeft,valLeft,BCRight,valRight,tlim,tord,physics,goal,x,h,k)
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
            end
        end
%     end
%     methods( Static=true)
        initializeexact(obj);
        computemoments(obj);
        computeprimalpseudo(obj);
        Z = unstructuredrecon(obj,u);
        er = reconplot(obj,Z);
        [uu,d] = updatesolution(obj,u);
        FI = computefluxintegral(obj,u);
        J=computefluxjacobian(obj,u);
        
    end
    
    methods(Access=private)
       
        poissoninitialize(obj);
        advectioninitialize(obj);
        burgersinitialize(obj);
    end
        
    
end

