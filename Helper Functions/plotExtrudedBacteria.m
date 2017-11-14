function plotExtrudedBacteria(time,cumExtrBact,extrSteps)

figure,scatter(time(extrSteps),cumExtrBact(extrSteps),400,'k.')
hold on
g = fit(time(extrSteps),cumExtrBact(extrSteps),'exp1');

smoothTime = linspace(0,time(max(extrSteps)));
plot(smoothTime,(g.a)*exp(g.b*smoothTime),'r','linewidth',2)
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}Number of Extruded Bacteria')
get(gca, 'XTick');
set(gca,'FontSize',20)
legend('\fontsize{20}Extruded Bacteria','Exp Fit','Location','northwest')
disp(['Rate of extrusion = ', num2str(g.b)])

end