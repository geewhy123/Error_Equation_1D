

V1 = 0.9558;
V2 = 0.3109;
V3 = 0.9418;
[U1,U2,U3]=toconservedvars(V1,V2,V3);


gam = 1.4;

[Uavg1 Uavg2 Uavg3] = computeroeavg(U1,U2,U3,U1,U2,U3);
AT = computeAtilde(Uavg1,Uavg2,Uavg3);
[rho,u,P] = toprimitivevars(Uavg1,Uavg2,Uavg3);
E = Uavg3/Uavg1; 
H = E+P/rho;

J = [ 0 1 0;
    0.5*(gam-3)*u^2 (3-gam)*u gam-1;
    u*(0.5*(gam-1)*u^2-H) H-(gam-1)*u^2 gam*u];

% assert(norm(J-AT)/norm(J) < 1e-10);




V4 = 0.9604;
V5 = 0.3155;
V6 = 0.9464;
[U4,U5,U6]=toconservedvars(V4,V5,V6);


[Uavg1 Uavg2 Uavg3] = computeroeavg(U1,U2,U3,U4,U5,U6);
AT = computeAtilde(Uavg1,Uavg2,Uavg3);
[F1 F2 F3] = conservedtoflux(U1,U2,U3);
[F4 F5 F6] = conservedtoflux(U4,U5,U6);
% assert(norm(AT*([U1 U2 U3]'-[U4 U5 U6]')- ([F1 F2 F3]'-[F4 F5 F6]'))/norm(AT*([U1 U2 U3]'-[U4 U5 U6]')) < 1e-10)


% assert(abs(det(AT))>1e-10)

% fprintf('Passed property 2,3,4')

AT
[U4 U1]
error('2')

0.5*([F1 F2 F3]'+[F4 F5 F6]')-0.5*AT*([U4; U5; U6]-[U1; U2; U3])
error('1')


