function [uu,d] = updatesoln(u,x,f,k,h,N,p,t,phys)

%UPDATESOLN Summary of this function goes here
%   Detailed explanation goes here


    
   switch t
       case 1
      [uu,d] = rk1(u,x,f,k,h,N,p,t,phys);

       case 2
           [uu,d] = rk2(u,x,f,k,h,N,p,t,phys);
           
       case 4
           [uu,d] = rk4(u,x,f,k,h,N,p,t,phys);
           
       case 5
           [uu,d] = rk5(u,x,f,k,h,N,p,t,phys);
           
       case 7
           [uu,d] = rk7(u,x,f,k,h,N,p,t,phys);
       otherwise
   end
      


end

