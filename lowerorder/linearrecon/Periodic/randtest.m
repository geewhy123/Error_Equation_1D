


close all 
clear all
Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
%Q = [2 2 3; 2 2 4 ;2 2 5];
% Q = [3 3 4];
% Q = [2 2 4; 2 2 6; 2 4 4;2 4 6; 2 6 4;2 6 6; 4 2 6; 4 4 6; 4 6 6 ]
% Q = [2 4 4; 2 6 4 ; 2 6 6 ]
%Q = [2 6 6];
% Q = [2 4 4; 2 6 4];
e20 = 0;
e40 =0;
% Q = [2 6 6]
[m,n] = size(Q);
for k= 1:m
for j = 1:1
% [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),1/3,100,10,7,'Poisson');
%[err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(80,Q(k,1),Q(k,2),Q(k,3),0,100,.3,7,'Advection');
 [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(80,Q(k,1),Q(k,2),Q(k,3),0,100,.2,7,'Burgers');
close all
e20 = e20+err2b(j);
% [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),1/3,100,10,7,'Poisson');
% [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),0,100,.3,7,'Advection');
 [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),0,100,.2,7,'Burgers');
close all
e40 = e40+err4b(j);
E(j) = log(err2b(j)/err4b(j))/log(2);
%E(j) = log(cverr2b(j)/cverr4b(j))/log(2);


end
Z(k) = mean(E);
end

Z
% E
% mean(E)
%log(e20/e40)/log(2)
save('Z.mat','Z')


% Q = [2 0 0 ; 3 0 0 ; 4 0 0 ; 5 0 0 ; 6 0 0 ];
% % Q = [6 0 0];
% e20 = 0;
% e40 =0;
% % Q = [2 6 6]
% [m,n] = size(Q);
% for k= 1:m
% for j = 1:5
% [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),1/3,100,.2,7,'Burgers');
% close all
% % e20 = e20+err2b(j);
% % [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),1/3,100,10,7,'Poisson');
% [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(60,Q(k,1),Q(k,2),Q(k,3),1/3,100,.2,7,'Burgers');
% close all
% % e40 = e40+err4b(j);
% %E(j) = log(err2b(j)/err4b(j))/log(2);
% E(j) = log(cverr2b(j)/cverr4b(j))/log(6/4);
% 
% 
% end
% Z(k) = mean(E);
% end

Z