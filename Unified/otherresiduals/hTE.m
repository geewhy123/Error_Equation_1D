function [ res ] = hTE( u,obj )
%HTE Summary of this function goes here
%   Detailed explanation goes here
x = obj.cellEndPoints;
N = length(x)-1;
xf = zeros(2*N+1,1)';
xf(1:2:end) = x;
xf(2:2:end-1) = (x(2:end)+x(1:end-1))/2;

 objf= obj.copy;
 objf.nCells= 2*N;
 objf.cellEndPoints = xf;
 objf.cellCentroids = [NaN (xf(2:end)+xf(1:end-1))/2 NaN];
 objf.cellWidths = [NaN xf(2:end)-xf(1:end-1) NaN];
 objf.computemoments();
 objf.initializeexact();
%  figure
% plot(objf.source,'*')
% error('2')
uf = zeros(2*N+2,1);
objf.computeprimalpseudo();
objf.primalJacobian = [];
% [Zf] = objf.unstructuredrecon(u,p,'solution')

 J = objf.computefluxjacobian(uf,'solution');
 J = J(2:end-1,2:end-1);
% ef = ones(2*N,1);
% dxf = dx/2;
% Af = spdiags([ef -2*ef ef],-1:1,2*N,2*N)/dxf^2;
% bcf = 0*Af;
% bcf(1,1) = -1/dxf^2;
% bcf(end,end) = -1/dxf^2;
% Af = Af+bcf;

% uf(1:2:end-1) = u;
% uf(2:2:end) = u;%(u(1:end-1)+u(2:end))/2;
% % uf(2:end-1) = (uf(1:end-2)+uf(3:end))/2;

% for i = 1:10
%     uf(1) = uf(2)/3;
%     for j = 2:length(uf)-1
%         uf(j) = (uf(j-1)+uf(j+1) )/2;
%     end
%     uf(end) = uf(end-1)/3;
% end
%try second order interp

uf = zeros(2*N,1);
u = u(2:end-1);
for i = 2:size(uf)-1
    uf(i) = 0.75*u(ceil(i/2));
    if mod(i,2) == 1
       uf(i) = uf(i) + 0.25*u(ceil(i/2)-1); 
    else
        uf(i) = uf(i) + 0.25*u(ceil(i/2)+1);
    end
end
uf(1) = u(1)/2;
uf(end) = u(end)/2;



% sf = (-4*pi^2/dxf)*(-1/(2*pi))*(cos(2*pi*xf(2:end))-cos(2*pi*xf(1:end-1)))';
sf = objf.source(2:end-1);

resf =J*uf-sf;
figure
plot(resf,'*')
 res = (resf(1:2:end-1)+resf(2:2:end))/2;
 res = [NaN; res;NaN];
  res = -res;
 
 
%  J*uf-sf
% res
%  error('1');

% %  res = A4*u-s;
% %  figure
% %  plot((x(1:end-1)+x(2:end))/2,resc,'o',(x(1:end-1)+x(2:end))/2,res,'-x')
% %  error('1')
% ee = A\-res;
% % max(abs(e_exact-ee))
% errerr = e_exact-ee;
% enorm = [mean(abs(errerr)) sqrt((1/N)*sum(errerr.^2)) max(abs(errerr))];
% disp(enorm);
% % ee


end

