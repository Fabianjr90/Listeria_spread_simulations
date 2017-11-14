function plotExtrusion(x,y,nlive,step,plotwidth,...
    hostCell,extruded)

scatter(x,y,10,'k','filled')
pbaspect([1 1 1]);
xlim([-plotwidth,plotwidth])
ylim([-plotwidth,plotwidth])
get(gca, 'XTick');
set(gca,'FontSize',20)
title(sprintf('Step %d, nlive: %d, extruded: %d',step, nlive,extruded ))
hold on
plot(hostCell(:,1),hostCell(:,2),'c--','LineWidth',2)
hold off
drawnow

end