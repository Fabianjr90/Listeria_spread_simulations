function plotConvexHull(data,areaConvHullTime,area_fit)

scatter(data(2:end,1),data(2:end,2),5,'ko','filled');
hold on;
u = plot(areaConvHullTime,area_fit,'r');
set(u,'linewidth',2)
xlim([0,max(areaConvHullTime)])
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}Area of Convex Hull')
get(gca, 'XTick');
set(gca,'FontSize',20)
legend('\fontsize{20}Raw Data','Model Fit','Location','northwest')

end