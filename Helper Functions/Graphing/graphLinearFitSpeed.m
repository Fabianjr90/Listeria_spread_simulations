function graphLinearFitSpeed(time)

xlim([1,max(time)])
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}Radial Distance')
get(gca, 'XTick');
set(gca,'FontSize',20)
legend('\fontsize{20}Raw data','Linear Fit','Location','northwest')

end