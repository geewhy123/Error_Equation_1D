[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Poisson');
assert(abs(errerr2-0.8347)<1e-2)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,1,7,'Advection');
assert(abs(errerr2-0.054)<1e-2)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.4,7,'Advection');
assert(abs(errerr2-0.0106)<1e-2)
close all
[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Burgers');
assert(abs(errerr2-0.0011)<1e-2)
close all
clc
printf('Passed All Tests')

