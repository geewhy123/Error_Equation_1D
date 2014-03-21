clearvars -except errerr10 errerr20 errerr40 errerr10b errerr20b errerr40b berrerr2 berrerr4 berrerr8 berrerr berrerr4 berrerr8 berrerrnon2 berrerrnon4
close all
A = [2 2 4;2 2 6;2 4 4;2 4 6; 2 6 4; 2 6 6; 4 2 6; 4 4 6; 4 6 6];
B = ['2 2 4';'2 2 6';'2 4 4';'2 4 6'; '2 6 4'; '2 6 6'; '4 2 6'; '4 4 6'; '4 6 6'];

[m,n] = size(A);
N = 80
figure
hold on
for j = 1:m
  % [errerr2b(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1,100,1,7,'Advection');
 
  %[errerr1b(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),1,100,10,7,'Poisson');
 %%% [errerr4(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Poisson');
 
  [berrerrnon8(j),x,cverr,exacterr(:,j),ee(:,j)] = errordriver(N,A(j,1),A(j,2),A(j,3),0,100,.1,7,'Burgers');
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

