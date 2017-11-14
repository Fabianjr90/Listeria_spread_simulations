function [area,xout,yout] = simConvexHullNEW(x,y,nlive,step,plotwidth)

[K,area] = convhull(x,y);
xout = x(K);
yout = y(K);

% scatplot(x,y,'circles');
% dscatter(x,y,'marker','o','msize',10);
% pbaspect([1 1 1]);
% xlim([-plotwidth,plotwidth])
% ylim([-plotwidth,plotwidth])
% get(gca, 'XTick');
% set(gca,'FontSize',20)
% title(sprintf('Step %d, nlive: %d',step, nlive))
% hold on
% u=plot(xout,yout,'r');
% set(u,'linewidth',2);
% hold off
% drawnow

end