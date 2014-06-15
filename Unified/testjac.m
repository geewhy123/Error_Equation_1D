U = [1 2 3];

gam = 1.4;


[rho,u,P] = toprimitivevars(U(1),U(2),U(3));
E = U(3)/U(1); 
H = E+P/rho;

J = [ 0 1 0;
    0.5*(gam-3)*u^2 (3-gam)*u gam-1;
    u*(0.5*(gam-1)*u^2-H) H-(gam-1)*u^2 gam*u];


U1 = U+ 1e-8*[1 0 0 ];
U2 = U+ 1e-8*[0 1 0 ];
U3 = U+ 1e-8*[0 0 1 ];


[F(1), F(2), F(3)]=conservedtoflux(U(1),U(2),U(3));

F

[F1(1), F1(2), F1(3)]=conservedtoflux(U1(1),U1(2),U1(3));

F1
[F2(1), F2(2), F2(3)]=conservedtoflux(U2(1),U2(2),U2(3));
[F3(1), F3(2), F3(3)]=conservedtoflux(U3(1),U3(2),U3(3));

dFdU(1:3,1) = (F1-F)'/1e-8;
dFdU(1:3,2) = (F2-F)'/1e-8;
dFdU(1:3,3) = (F3-F)'/1e-8

J


J*U'