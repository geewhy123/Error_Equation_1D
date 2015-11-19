function [uu,d] = updateerror(obj,u,TT,j)%eqn,u,x,f,k,h,N,p,t,phys,time,Rsp)

%UPDATESOLN Summary of this function goes here
%   Detailed explanation goes here

x = obj.cellCentroids;
f = obj.residual(:,j)*-1;
k = obj.tStep;
h = obj.cellWidths;
N = obj.nCells;
p = obj.qOrder;
phys = obj.physics;    
t = obj.tOrder;
% time = NaN;
time = TT;
% Rsp = NaN;
Rsp = obj.Rsp;
BCLeft = obj.bcLeftType;
uL = obj.bcLeftVal;
BCRight = obj.bcRightType;
uR = obj.bcRightVal;
   switch t
       case 1
       [uu,d] = rk1('error',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       case 2
%            [uu,d] = rk2(u,x,f,k,h,N,p,t,phys);
           [uu,d] = rk2('error',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       case 3
         [uu,d] = rk3('error',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
%              [uu,d] = rk4('error',u,x,f,k,h,N,p,phys,time,Rsp,'D',1,'D',1,obj);
       case 4
          
%            [uu,d] = rk4(u,x,f,k,h,N,p,t,phys);
 [uu,d] = rk4('error',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
%              [uu,d] = rk4('error',u,x,f,k,h,N,p,phys,time,Rsp,'D',1,'D',1,obj);
       case 5
           [uu,d] = rk5(u,x,f,k,h,N,p,t,phys);
           
       case 7
           
%          [uu,d] = rk7('error',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
           [uu,d] = rk7('error',u,x,f,k,h,N,p,phys,time,Rsp,'D',1,'D',1,obj);
       otherwise
   end
      


end

