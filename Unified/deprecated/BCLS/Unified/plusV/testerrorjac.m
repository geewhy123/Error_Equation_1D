syms u1 u2 u3 e1 e2 e3 gam


g = [e2; 
(gam-1)*e3+((3-gam)/2)*((e2+u2)^2/(e1+u1) - u2^2/u1) ;  
gam*((e2+u2)*(e3+u3)/(e1+u1) -u2*u3/u1) - ((gam-1)/2)*((e2+u2)^3/(e1+u1)^2-u2^3/u1^2) ];

J = [diff(g(1),e1) diff(g(1),e2) diff(g(1),e3) ; diff(g(2),e1) diff(g(2),e2) diff(g(2),e3); diff(g(3),e1) diff(g(3),e2) diff(g(3),e3) ]
    
syms x1 x2 x3

A = [0 1 0; ((gam-3)/2)*(x2^2/x1^2) (3-gam)*(x2/x1) gam-1; -gam*x3*x2/x1^2+(gam-1)*x2^3/x1^3 x3*gam/x1+(1-gam)*(3/2)*x2^2/x1^2 x2*gam/x1]


