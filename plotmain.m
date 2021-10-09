%%
close all;

clear all;
clc

%%
load('consenuswithoutcommunication');
figure(1)
plot([0:0.01:4.01],V,'linewidth',2)
xlim([0 4])
set(gca,'XTick',[0:1:4]);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlabel('Time (seconds)')
ylabel('DGU voltage')
hold on
plot([0:0.01:4.01], [ones(201,1)*Vmax1;ones(201,1)*Vmax2],'b--','linewidth',2);
hold on
plot([0:0.01:4.01], [ones(201,1)*Vmin1;ones(201,1)*Vmin2],'b--','linewidth',2);
grid on
legend('$V_1$','$V_2$','$V_3$','$V_4$','$V_5$','$V_6$','Location','north','Interpreter','latex','Orientation','horizontal')
hold off plots
%%
figure(2);
plot([0:0.01:4.01],Ifilt*D,'linewidth',2);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
leg=legend('$\frac{I_{t1}}{\bar{I}_{t1}}$','$\frac{I_{t2}}{\bar{I}_{t2}}$','$\frac{I_{t3}}{\bar{I}_{t3}}$','$\frac{I_{t4}}{\bar{I}_{t4}}$','$\frac{I_{t5}}{\bar{I}_{t5}}$','$\frac{I_{t6}}{\bar{I}_{t6}}$','Location','north','Orientation','horizontal','Interpreter','latex');
leg.FontSize = 20;
xlabel('Time (seconds)')
xlim([0 4]);
set(gca,'XTick',[0:1:4]);
ylim([4,6]);
ylabel('Weighted filter currents')
grid on
hold off plots

%%
figure(3);
plot([0:0.01:4.01],V*D^(-1)*ones(6,1),'r','linewidth',2);
hold on
plot([0:0.01:4.01], (ones(1,6)*D^(-1)*Vref)*ones(402,1),'b--','linewidth',2.5);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman');
xlim([0 4]);
set(gca,'XTick',[0:1:4]);
xlabel('Time (seconds)')
ylabel('Voltage regulation')
grid on
hold off plots