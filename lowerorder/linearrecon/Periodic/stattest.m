% % % 
% % % 
% % % 
% % % close all 
% % % clear all
% % % Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]
% % % %Q = [2 2 3; 2 2 4 ;2 2 5];
% % % % Q = [3 3 4];
% % % % Q = [2 2 4; 2 2 6; 2 4 4;2 4 6; 2 6 4;2 6 6; 4 2 6; 4 4 6; 4 6 6 ]
% % % % Q = [2 4 4; 2 6 4 ; 2 6 6 ]
% % % %Q = [2 6 6];
% % % % Q = [2 4 4; 2 6 4];
% % % % Q = [3 6 5;]
% % % e20 = 0;
% % % e40 =0;
% % % % Q = [2 6 6]
% % % % Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6;]
% % % % Q = [2 2 4;]
% % % % Q = [2 2 3; 2 2 4];
% % % 
% % % [m,n] = size(Q);
% % % for k= 1:m
% % % for j = 1:500
% % % % [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),1/3,100,10,7,'Poisson');
% % % [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),1/3,100,.3,7,'Advection');
% % % %  [err2b(j),x,cverr2b(j),exacterr2b(:,j),ee2b(:,j)] = errordriver(20,Q(k,1),Q(k,2),Q(k,3),0,100,.2,7,'Burgers');
% % % close all
% % % e20 = e20+err2b(j);
% % % % [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),1/3,100,10,7,'Poisson');
% % %  [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),1/3,100,.3,7,'Advection');
% % % %  [err4b(j),x,cverr4b(j),exacterr4b(:,j),ee4b(:,j)] = errordriver(40,Q(k,1),Q(k,2),Q(k,3),0,100,.2,7,'Burgers');
% % % close all
% % % e40 = e40+err4b(j);
% % % E(j,k) = log(err2b(j)/err4b(j))/log(2);
% % % %E(j) = log(cverr2b(j)/cverr4b(j))/log(2);
% % % 
% % % 
% % % % normE = (E-mean(E))./sqrt(var(E))
% % % end
% % %  Z(k) = mean(E(:,k));
% % % end
% % % 
% % % Z
% % % % E
% % % % mean(E)
% % % %log(e20/e40)/log(2)
% % % save('Z.mat','Z')
% % % 
% % % hist(E)
close all
  Q = [2 2 3; 2 2 4; 2 2 5; 2 2 6; 2 3 3; 2 3 4; 2 3 5; 2 3 6; 2 4 3; 2 4 4; 2 4 5; 2 4 6; 2 5 3; 2 5 4; 2 5 5; 2 5 6; 2 6 3; 2 6 4; 2 6 5; 2 6 6; 3 2 4; 3 2 5; 3 2 6; 3 3 4; 3 3 5; 3 3 6; 3 4 4; 3 4 5; 3 4 6; 3 5 4; 3 5 5; 3 5 6; 3 6 4; 3 6 5; 3 6 6; 4 2 5; 4 2 6; 4 3 5; 4 3 6; 4 4 5; 4 4 6; 4 5 5 ; 4 5 6; 4 6 5; 4 6 6; 5 2 6; 5 3 6 ;5 4 6 ; 5 5 6; 5 6 6  ]

p = Q(:,1);
q = Q(:,2);
r = Q(:,3);

fit = min(p+q,min(r,p-1)+1)+max((r-p).*(q==r),(p==r-1));

plot(fit,'*-')
hold on
boxplot(E)


 q = [223; 224; 225; 226; 233; 234; 235; 236; 243; 244; 245; 246; 253; 254; 255; 256; 263; 264; 265; 266; 324; 325; 326; 334; 335; 336; 344; 345; 346; 354; 355; 356; 364; 365; 366; 425; 426; 435; 436; 445; 446; 455 ; 456; 465; 466; 526; 536 ;546 ;  556; 566  ];

text = strtrim(cellstr(num2str(q)));
% % xlabel(text)
 set(gca,'XTick',1:1:50,'XTickLabel',text)
ylabel('average(100) order of error convergence')

figure


for j = 1:50
    X(j) = chi2gof(E(:,j));
%    if(chi2gof(E(:,j))==1)
%       error('1') 
%    end


subplot(10,5,j)
% hist(E(:,j))
if(X(j)==1)
set(get(gca,'child'),'FaceColor','r');
end

 xlabel(num2str(Q(j,:)))
% if (k==1)
%    legend('exact error','computed error'); 
% end
hold on
% subplot(2,1,k)
% plot(x,exacterr(:,k),'o',x,ee(:,2),'*');
% figure
% e = E(:,j);
% mu = mean(e);
% s = sqrt(var(e));
% x = mu-5*s:0.0001:mu+5*s;
% y = normpdf(x,mu,s);
% plot(x,y);
% error('1');
histfit(E(:,j))

% t(j)=ttest(E(:,j),fit(j));
x = [fit(j),fit(j)];
y = [0,100];
plot(x,y,'*b-')




e = E(:,j);
e = (e-mean(e))/sqrt(var(e));
[h(j),p(j)]=kstest(e);

end
% figure
% hist(E(:,1))
% set(get(gca,'child'),'FaceColor','r');






fprintf('confidence 0.95')
 save('test.mat','E')
