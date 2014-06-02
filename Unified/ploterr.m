clearvars -except errerr10 errerr20 errerr40 errerr80 errerr10b errerr20b errerr40b berrerr2 berrerr4 berrerr8 berrerr berrerr4 berrerr8 berrerrnon2 berrerrnon4 berrerrlin2 berrerrlin4 err2b cverr2b err4b
close all
A = [2 2 4;2 2 6;2 4 4;2 4 6; 2 6 4; 2 6 6; 4 2 6; 4 4 6; 4 6 6];
B = ['2 2 4';'2 2 6';'2 4 4';'2 4 6'; '2 6 4'; '2 6 6'; '4 2 6'; '4 4 6'; '4 6 6'];

A = [2 2 5; 2 3 4; 2 3 5; 3 3 5; 3 4 5; 3 5 5 ; 4 2 5; 4 4 5; 4 5 6];
A = [3 2 4; 3 2 5; 3 2 6; 3 3 4; 4 3 5; 4 3 6; 4 2 5; 4 2 6; 5 2 6];
 A = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5 ; 2 3 6; 2 4 3;2 4 4; 2 4 5 ;2 4 6 ;2 5 3 ;2 5 4 ;2 5 5 ; 2 5 6;2 6 3; 2 6 4; 2 6 5; 2 6 6;]%3 2 4; 3 2 5 ; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4 ; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 

% % % %  B = ['2 2 3'; '2 2 4'; '2 2 5';'2 2 6'; '2 3 3'; '2 3 4'; '2 3 5' ; '2 3 6'; '2 4 3';'2 4 4'; '2 4 5' ;'2 4 6' ;'2 5 3' ;'2 5 4' ;'2 5 5' ; '2 5 6';'2 6 3'; '2 6 4'; '2 6 5'; '2 6 6';]%3 2 4; 3 2 5 ; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4 ; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 
A = [2 2 4;2 2 6;2 4 4;2 4 6; 2 6 4; 2 6 6; 4 2 6; 4 4 6; 4 6 6];
B = ['2 2 4';'2 2 6';'2 4 4';'2 4 6';  '4 2 6'; '4 4 6'; ];
A = [  2 2 4;  2 2 6;2 4 4;2 4 6; 4 2 6; 4 4 6; ];
% A = [2 2 4;2 4 4];
[m,n] = size(A);

N = 20;
figure
hold on
for j = 1:m
%    [err4b(j),x,cverr4b(j),exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1/10,100,.1,7,'Advection');
 
%   [errerr80(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1/3,'D',0,'D',0,10,7,'Poisson','SS');
  [errerr20(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,'D',1-tanh(0),'D',1-tanh(1/2),10,7,'BurgersMod','SS');
  
  
 %%% [errerr4(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Poisson');
 
%  [berrerrlin8(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Burgers');
close all

end


for k = 1:m
subplot(3,3,k)
plot(x,exacterr(:,k),'o',x,ee(:,k),'*');

xlabel(B(k,:))
if (k==1)
   legend('exact error','computed error'); 
end
hold on
% subplot(2,1,k)
% plot(x,exacterr(:,k),'o',x,ee(:,2),'*');
end

