[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Poisson');
assert(abs(errerr2-0.8347)<1e-3)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,1,7,'Advection');
assert(abs(errerr2-0.054)<1e-3)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.4,7,'Advection');
assert(abs(errerr2-0.0106)<1e-3)
close all
%[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Burgers');
%assert(abs(errerr2-0.0011)<1e-3)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,4,4,1/10,100,10,7,'Poisson');
assert(abs(errerr2-0.0029)<1e-3)
close all

clc
printf('Passed All Tests')



% Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
% 
% [m,n] = size(Q);
% j=1;
% for k= 1:m
% 
% 
% [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),0,100,.1,7,'Advection');
% close all
% 
% 
% [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),0,100,.1,7,'Advection');
% close all
% 
% E(k) = log(err2b(j)/err4b(j))/log(2);
% 
% 
% end
% E