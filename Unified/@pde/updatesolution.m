function [uu,d] = updatesolution(obj,u)%eqn,u,x,f,k,h,N,p,t,phys,time,Rsp)

%UPDATESOLN Summary of this function goes here
%   Detailed explanation goes here

x = obj.cellCentroids;
f = obj.source;
k = obj.tStep;
h = obj.cellWidths;
N = obj.nCells;
p = obj.pOrder;
phys = obj.physics;    
t = obj.tOrder;
time = NaN;
Rsp = NaN;
BCLeft = obj.bcLeftType;
uL = obj.bcLeftVal;
BCRight = obj.bcRightType;
uR = obj.bcRightVal;
   switch t
       case 'rk1'
           %       [uu,d] = rk1(u,x,f,k,h,N,p,t,phys);
           [uu,d]= rk1('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       case 'rk2'
           %            [uu,d] = rk2(u,x,f,k,h,N,p,t,phys);
           [uu,d] = rk2('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
        
       case 'rk3'
           [uu,d] = rk3('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       case 'rk4'
          
           [uu,d] = rk4('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
           
       case 'rk5'
           [uu,d] = rk5(u,x,f,k,h,N,p,t,phys);
           
       case 'rk7'
           if(strcmp(obj.physics,'LinearSystem')~=1)
           [uu,d] = rk7('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
           else
           [uu,d] = rk7LinearSystem('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
           end
       case 'irk1'
           [uu,d] = irk1('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       
       case 'irk2'
           [uu,d] = irk2('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
          case 'irk4'
           [uu,d] = irk4('solution',u,x,f,k,h,N,p,phys,time,Rsp,BCLeft,uL,BCRight,uR,obj);
       
       otherwise
           error('1')
   end
      


end

