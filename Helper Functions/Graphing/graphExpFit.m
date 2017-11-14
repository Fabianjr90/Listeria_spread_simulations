function graphExpFit(time,bacteria_number,g,step)
plot(time(1:(step-1)),(g.a)*exp(g.b*time(1:(step-1))),'r','linewidth',2)
hold on;
scatter(time(1:10:(step-1)),bacteria_number(1:10:(step-1)),200,'k.')
hold off;
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}Number of Bacteria')
get(gca, 'XTick');
set(gca,'FontSize',20)
legend('\fontsize{20}Exponential Fit','nlive','Location','northwest')

end