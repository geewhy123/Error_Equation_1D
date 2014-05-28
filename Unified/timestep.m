function [ phi]  = timestep( equation, Z,f,k,h,N,p,phys,time,Rsp,Zu,val,obj)
%   Detailed explanation goes here

switch equation
    case 'solution'
        phi=reconfluxsoln( Z,f,h,N,p,phys,time);
    case 'error'
        obj.errorSource = -1*val;
         phi = obj.computefluxintegral(Z,'error');

%           phi=reconfluxerror( Z,f,k,h,N,p,phys,time,Rsp,Zu,val,obj);
%           phi
%           error('1')
    otherwise
        error('2')
    

end

