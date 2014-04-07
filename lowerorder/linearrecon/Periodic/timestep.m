function [ output_args ] = timestep( equation)
%TIMESTEP Summary of this function goes here
%   Detailed explanation goes here


switch equation
    case 'solution'
        reconfluxsoln
    case 'error'
        reconfluxerror
    otherwise
        error('2')
    

end

