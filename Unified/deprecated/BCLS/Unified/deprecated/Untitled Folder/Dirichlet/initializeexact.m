function [ u0,ue,f] = initializeexact(physics,N,x,h,tlim )
%INITIALIZEEXACT Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(physics,'Burgers')==1)
   [u0,ue,f]=burgersinitialize(N,x,h,tlim);

elseif(strcmp(physics,'Advection')==1)
     [u0,ue,f]=advectioninitialize(N,x,h,tlim);
elseif(strcmp(physics,'Poisson')==1)
     [u0,ue,f]=poissoninitialize(N,x,h,tlim);
else
    error('1')
end


end

