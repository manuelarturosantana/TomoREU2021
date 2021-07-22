figure(1)

load LSQRcontrola.mat
hold on
plot(iterInfo.xErrors,'-d','LineWidth',1.5)

load IrnControl.mat
plot(iterInfo.xErrors,'-^','LineWidth',1.5)

load FistaControl.mat
plot(iterInfo.xErrors,'-*','LineWidth',1.5)

legend('LSQR','IRN','FISTA')
xlabel('Iterations','FontSize',20)
ylabel('Relative Error','FontSize',20)
ax = gca;
ax.FontSize = 16;
