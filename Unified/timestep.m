function [ phi]  = timestep( equation, Z,f,k,h,N,p,phys,time,Rsp,Zu,val,obj)
%   Detailed explanation goes here

switch equation
    case 'solution'
        phi=reconfluxsoln( Z,f,h,N,p,phys,time);
    case 'error'
%         phi = obj.computefluxintegral(Z,
         phi=reconfluxerror( Z,f,k,h,N,p,phys,time,Rsp,Zu,val,obj);
    otherwise
        error('2')
    

end

