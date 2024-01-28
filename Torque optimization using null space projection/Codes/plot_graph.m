o=1:1:628;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(o,opti,'LineWidth',1.5);
title("Optimization",'Interpreter','latex');
xlabel('Sample Numbers','Interpreter','latex');
ylabel('$\sqrt{\tau^\prime \tau}$', 'Interpreter', 'latex', 'FontSize', 14);
ylim([0 2]);
set(gca,'FontSize',18);
grid minor;
hold on
plot(o,eff,'LineWidth',1.5);
legend('With Optimisation','Without Optimisation','Interpreter','latex');
hold off