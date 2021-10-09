clear all

%% Network definition
com=1; %com=0 means no communication
N=6;
M=7;
% Incidence Matrix
B =    [1         1    0    0    -1         0        0 ;
-1        0    -1   0     0         0         0 ;
0         -1   0    -1    0         0         0 ;
0         0    1    1     0         -1        0 ;
0         0    0    0     0         1         1 ;
0         0    0    0     1         0        -1 ];


% Parameters of the filter

C_t=diag([2.2,1.9,1.7,2.5,2.0,3.0])*1e-3;
R_t=diag([0.2,0.3,0.1,0.5,0.4,0.6]);
L_t=diag([1.8,2.0,2.2,3.0,1.2,2.5])*1e-3;

%Parameters of the Network
R=diag([05,07,04,06,10,08,08])*0.1;
L=diag([2.1,1.8,2.3,1.0,1.8,2.5,2.0])*1e-4;

%load parameters
Y_l=diag([20,20,40,20,30,40])^-1;
%Y_l=diag([10,10,10,10,10,10])^-1;
%P_l=diag([100,100,50,100,70,50]);
P_l=100*eye(N);
I_l=[2;4.5;2.5;3.5;2.75;1];
expo=0;

%communication Laplacian
% Lc=10*[2 -1 0 0 0 -1;
% -1 3 -1 -1 0 0;
% 0 -1 2 0 -1 0;
% 0 -1 0 1 0 0;
% 0 0 -1 0 2 -1;
% -1 0 0 0 -1 2];

Lc=10*[1 0 0 0 -1 0;
0 2 -1 -1 0 0;
0 -1 2 0 -1 0;
0 -1 0 1 0 0;
-1 0 -1 0 3 -1;
0 0 0 0 -1 1];
Le=B*R^(-1)*B';
%D=eye(N);
D=diag([1/1.5,1/1.08,1/1.2,1/1.15,1/1,1/1.15]);


%Controller Gains 
k1=-eye(6);
k2=R_t/2;
k3=L_t^(-1)*R_t/2;

alpha=L_t^(-1)*(k1-eye(N));
beta=L_t^(-1)*(k2-R_t);
gamma=L_t^(-1)*k3;


%%
Vref=[50.1,50.5,51,51.5,49.4,50.4]';
Ist=inv(D);
Lt=Ist-sum(sum(Ist))^(-1)*Ist*ones(N,N)*Ist;
Itilde=[-Lt*(Ist)^(-1)*I_l; ones(1,N)*Ist*Vref];
Ltilde=[Le+Lt*Ist^(-1)*Y_l;ones(1,N)*Ist];
Lttilde=[Lt*Ist^(-1);zeros(1,N)];
Vstar=pinv(Ltilde)*Itilde;
criticalpower=(diag(Vstar)^-1)*4*pinv(Ltilde)*Lttilde*diag(Vstar)^-1;
Delta=norm(criticalpower*sum(P_l)',inf);
deltaminus=(1-sqrt(1-Delta))/2;
deltaplus=(1+sqrt(1-Delta))/2;
upper=(1+deltaminus)*Vstar;
lower=(1-deltaminus)*Vstar;
Vmax=max(upper);
Vmin=min(lower);

%% Simulation initiliazation

t_init = 5;      % seconds for initilization time

t_sim = 15;           % seconds per simulation
T = 7;                % number of simulation iteration
Tot_sim = t_sim*T;    % seconds of total simulation time

%%%%%%%%%%%%%%%%%%%%%%%% Initialization
tr=[0:0.01:t_init];
display('*** Simulation started ***');

if com==1
%inval=[50.6922098927866,50.2075838536882,50.7350519471104,50.1197929294029,50.2490159836503,50.8205686894574,7.30064967624129,5.25646784787322,5.84051993096869,5.59716479723273,4.86709989081456,5.59716501210743,1.81674116955864,1.35717960239422,4.50084914376505,1.25279149450424,0.598640519867463,0.853992912352441,-0.266903533363419,0.234331217464776,0.0440265402093501,0.393054967648322,-0.178071945688614,-0.226437246270414,0.969254136137057,-0.0611959119775750,-0.219475812936127,-1.02543878924287,0.128371423372203,0.161526880199520,-0.714438546676929];
inval=[50.0000000050207,50.0000000011690,49.9999999974052,50.0000000039733,49.9999999999966,49.9999999937407,8.50000000011511,9.42857142859462,6.24999999995382,9.00000000005315,7.16666666666254,7.24999999989962,1.81530000000011,1.35219047619050,4.41374999999996,1.22700000000008,0.608599999999993,0.851458333333210,zeros(1,N),-7.16465681730959e-08,-4.02189394473282e-07,-3.80823818386207e-08,6.20127760619419e-07,-1.01942996012030e-06,9.13692791042091e-08,-9.16987288958291e-08];
else
inval=[50.8294801928396,50.4543171779561,51.0290091250504,50.3409531591373,50.0700755738917,50.0375077199018,7.97915448933590,5.74499135347782,6.38332387575893,6.11735184311639,5.31943648582484,6.11735216472144,1.81796247829244,1.35815664942119,4.50204331241937,1.25435205569774,0.599183323780726,0.855293380145715,0.750329109955723,-0.285030821717194,-0.283407805312421,-1.14677056123640,-0.791953571732207,-0.338599881071088,0.0407133124102873];    
end
networksimulationnahata
%% Plotting results
close all;
figure(1);
plot(t,y(:,1:N),'linewidth',2);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
hold on
if expo==0 && det(P_l)~=0 && com==1
plot(t, ones(size(t))*Vmax,'b--','linewidth',2);
hold on
plot(t, ones(size(t))*Vmin,'b--','linewidth',2);
end
xlabel('Time (seconds)')
ylabel('DGU voltage')
grid on
legend('$V_1$','$V_2$','$V_3$','$V_4$','$V_5$','$V_6$','Location','north','Interpreter','latex','Orientation','horizontal')
hold off plots

figure(2);
plot(t,y(:,N+1:2*N),'linewidth',2);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
legend('$I_{t1}$','$I_{t2}$','$I_{t3}$','$I_{t4}$','$I_{t5}$','$I_{t6}$','Location','southeast','Orientation','horizontal','Interpreter','latex')
xlabel('Time (seconds)')
ylabel('Filter currents')
grid on
hold off plots

if com==1
figure(3);
plot(t,y(:,N+1:2*N)*D,'linewidth',2);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
leg=legend('$\frac{I_{t1}}{\bar{I}_{t1}}$','$\frac{I_{t2}}{\bar{I}_{t2}}$','$\frac{I_{t3}}{\bar{I}_{t3}}$','$\frac{I_{t4}}{\bar{I}_{t4}}$','$\frac{I_{t5}}{\bar{I}_{t5}}$','$\frac{I_{t6}}{\bar{I}_{t6}}$','Location','southeast','Orientation','horizontal','Interpreter','latex');
leg.FontSize = 20;
xlabel('Time (seconds)')
ylabel('Weighted filter currents')
grid on
hold off plots

figure(4);

plot(t,y(:,1:N)*D^(-1)*ones(N,1),'r','linewidth',2);
hold on
plot(t, (ones(1,N)*D^(-1)*Vref)*ones(size(t)),'b--','linewidth',2.5);
set(gcf, 'Position',[189, 611,560,310]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlabel('Time (seconds)')
ylabel('Voltage regulation')
grid on
hold off plots
end
% figure(3);
% plot(t,y(:,19:25),'linewidth',1.5);
% set(gcf, 'Position',[189, 611,560,310]);
% set(gca, 'FontName', 'Times New Roman')
% legend('I_{1}','I_{2}','I_{3}','I_{4}','I_{5}','I_{6}','I_{7}','Location','southeast')
% xlabel('Time (seconds)')
% ylabel('Power (KW)')
% grid on
% hold off plots

% Reference Voltage
% Vref=50*ones(N,1)-0.2*rand(N,1);
% Ldash=inv(sum(sum(inv(D))))*inv(D)*ones(N,N)-eye(N);
% Adash=[Lp-Ldash*Y_l;ones(1,N)*D^-1];
% Istar=[Ldash*I_l; ones(1,N)*D^-1*Vref];
% Vstar1=pinv(Adash)*Istar;
% criticalpower1=diag(Vstar1)^-1*pinv(Adash)*[-Ldash;zeros(1,N)]*diag(Vstar1)^-1;
% Delta1=norm(criticalpower1*sum(P_l)',inf);
% deltaminus=(1-sqrt(1-Delta1))/2;
% deltaplus=(1+sqrt(1-Delta1))/2;
% upper1=(1+deltaminus)*Vstar1;
% lower1=(1-deltaminus)*Vstar1;
