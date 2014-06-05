
tic


[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,2,4,0,'P',0,'P',0,10,7,'Poisson','SS');
assert(abs(errerr2-0.014130429330857)/errerr2 < 0.001)
% [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,2,4,0,'P',0,'P',0,10,7,'Poisson','SS');
% assert(abs(errerr2-6.829423497046346e-04)/errerr2 < 0.001) 


[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,4,4,1/3,'P',0,'P',0,10,7,'Poisson','SS');
assert(abs(errerr2-0.002502354702354)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,2,4,1/3,'D',0,'D',0,10,7,'Poisson','SS');
assert(abs(errerr2-0.304897658837351)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,4,4,0,'D',0,'D',0,10,7,'Poisson','SS');
assert(abs(errerr2-0.001964073234010)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,2,4,0,'P',0,'P',0,0.3,7,'Advection','TimeAccurate');
assert(abs(errerr2- 0.010395717998870)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,4,4,1/3,'P',0,'P',0,0.3,7,'Advection','TimeAccurate');
assert(abs(errerr2-0.006105285299128)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,6,6,1/3,'P',0,'P',0,10,7,'Poisson','SS');
assert(abs(errerr2-3.915450275221799e-04)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,4,6,1/3,'D',0,'D',0,10,7,'Poisson','SS');
assert(abs(errerr2-0.064614874341152)/errerr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,4,0,0,1/3,'D',1-tanh(0),'D',1-tanh(1/2),10,7,'BurgersMod','SS');
assert(abs(cverr2-5.683313096888770e-06)/cverr2 < 0.001)

[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,4,4,1/3,'D',1-tanh(0),'D',1-tanh(1/2),10,7,'BurgersMod','SS');
% assert(abs(errerr2-5.774687156486830e-06)/errerr2 < 0.001)
assert(abs(errerr2-5.683313774878044e-06)/errerr2 < 0.001)


% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,4,1/3,'D',0,'D',0,0.3,7,'Advection','TimeAccurate');
% assert(abs(errerr2-0.304897658837351)/errerr2 < 0.001)

% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,4,4,0,'D',0,'D',0,0.3,7,'Advection','TimeAccurate');
% assert(abs(errerr2-0.001964073234010)/errerr2 < 0.001)
[errerr2,x,cverr2,exacterr,ee,te]=errordriver(10,2,0,0,0,'F',0,'D',1,0.5,7,'Advection','TimeAccurate');
assert(abs(cverr2-0.662135389183426)/cverr2 <0.001)

fprintf('Passed no 6th order')

% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Poisson');
% % % assert(abs(errerr2-0.8347)<1e-4)
% % % close all
%     [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,2,6,0,'P',0,'P',0,.1,7,'Poisson','TimeAccurate');
% assert(abs(errerr2-0.04129)<1e-5)
% close all
% 
% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,1,7,'Advection');
% assert(abs(errerr2-0.054)<1e-3)
% close all
% % % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.4,7,'Advection');
% % % assert(abs(errerr2-0.0106)<1e-4)
% % % close all
% [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,2,6,0,100,.6,7,'Advection');
% assert(abs(errerr2-7.0808e-4)<1e-8)
% close all
% 
% %[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Burgers');
% %assert(abs(errerr2-0.0011)<1e-3)
% close all
% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,4,4,1/10,100,10,7,'Poisson','SS');
% assert(abs(errerr2-0.0029)<1e-4)
% close all
% 
% 
% [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,4,1/10,100,10,7,'Poisson','FI');
% assert(abs(errerr2)<1e-12)


toc



% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % count = 0;
% % % while(1)
% % %     [errerr2,x,cverr2,exacterr,ee]=errordriver(40,2,6,6,1/3,100,10,7,'Poisson');
% % %     if(abs(errerr2) > 1e-6) 
% % % %    [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,6,6,1/3,100,10,7,'Poisson');
% % % %       if(abs(errerr2) > 2.5e-5) 
% % %         count
% % %         error('1')
% % %     end
% % %     count = count +1;
% % % end
% % 
% % 
% % 
% % tic
% % % % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Poisson');
% % % % assert(abs(errerr2-0.8347)<1e-4)
% % % % close all
% %     [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,2,6,0,'P',0,'P',0,.1,7,'Poisson','TimeAccurate');
% % assert(abs(errerr2-0.04129)<1e-5)
% % close all
% % 
% % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,1,7,'Advection');
% % assert(abs(errerr2-0.054)<1e-3)
% % close all
% % % % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.4,7,'Advection');
% % % % assert(abs(errerr2-0.0106)<1e-4)
% % % % close all
% % [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,2,6,0,100,.6,7,'Advection');
% % assert(abs(errerr2-7.0808e-4)<1e-8)
% % close all
% % 
% % %[errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,6,0,100,.1,7,'Burgers');
% % %assert(abs(errerr2-0.0011)<1e-3)
% % close all
% % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,4,4,1/10,100,10,7,'Poisson','SS');
% % assert(abs(errerr2-0.0029)<1e-4)
% % close all
% % 
% % 
% % [errerr2,x,cverr2,exacterr,ee]=errordriver(10,2,2,4,1/10,100,10,7,'Poisson','FI');
% % assert(abs(errerr2)<1e-12)
% % close all
% % 
% % 
% % clc
% % fprintf('Passed All Tests')
% % 
% % 
% % % % [errerr2,x,cverr2,exacterr,ee]=errordriver(20,2,4,6,0,100,.1,7,'Burgers');
% % % % assert(abs(errerr2-2.1328e-7)<1e-7)
% % % % close all
% % 
% % % % % [errerr2,x,cverr2,exacterr,ee]=errordriver(40,2,4,6,0,100,.2,7,'Burgers');
% % % % % assert(abs(errerr2-2.1360e-8)<1e-12)
% % % % % close all
% % % % % 
% % % % % fprintf('Passed Burgers')
% % 
% % % Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
% % % 
% % % [m,n] = size(Q);
% % % j=1;
% % % for k= 1:m
% % % 
% % % 
% % % [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),0,100,.1,7,'Advection');
% % % close all
% % % 
% % % 
% % % [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),0,100,.1,7,'Advection');
% % % close all
% % % 
% % % E(k) = log(err2b(j)/err4b(j))/log(2);
% % % 
% % % 
% % % end
% % % E
% % 
% % 
% % toc
