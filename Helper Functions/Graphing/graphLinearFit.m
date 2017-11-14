function graphLinearFit(time)

xlim([1,max(time)])
xlabel('\fontsize{25}Time')
ylabel('\fontsize{25}< r^{2} >')
get(gca, 'XTick');
set(gca,'FontSize',20)
legend('\fontsize{20}Raw data','Sanity Check','Location','northwest')

end