clear
figure(1)

xvalues = 5:20;
load bNoiseControl20.mat
hold on
plot(xvalues,xErrors(5:end),'-*','LineWidth',2)

load CScontrol.mat
plot(xvalues,xErrors(5:20),'-d','LineWidth',2)

load AAcontrol.mat %color 'k' etc
plot(xvalues,xErrors(5:20),'-^','LineWidth',2)

load ITcontrolb.mat
plot(xvalues,xErrors(5:20),'-o','LineWidth',2)

legend('BCD','CS','AA','IT')
xlabel('Iterations','FontSize',20)
ylabel('Relative Error','FontSize',20)
title('Image Error','FontSize',16)
ax = gca;
ax.FontSize = 16;