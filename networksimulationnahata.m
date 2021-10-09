% new filter parameters 

net_par = struct;
net_par.P_l  = P_l;
net_par.Vref= Vref;
net_par.C_t = C_t;
net_par.I_l = I_l;
net_par.M = M;
net_par.N = N;

if com==1
X=[-C_t^(-1)*Y_l,C_t^(-1),zeros(N,N),zeros(N,N),-C_t^(-1)*B;
    alpha,beta,gamma,alpha*D*Lc, zeros(N,M);
    -eye(N),zeros(N),zeros(N),-D*Lc,zeros(N,M);
    zeros(N,N),Lc*D, zeros(N,N),zeros(N,N),zeros(N,M);
   L^(-1)*B',zeros(M,N),zeros(M,N),zeros(M,N),-L^(-1)*R];


[t,y]=ode45(@(t,y) ZIPloads_par(t,y,X,net_par,expo),tr,inval);

% X=[-C_t^(-1)*Y_l,C_t^(-1),zeros(N,N),zeros(N,N),-C_t^(-1)*B;
%     alpha,beta,gamma,alpha*D*Lc*0, zeros(N,M);
%     -eye(N),zeros(N),zeros(N),eye(N),zeros(N,M);
%     zeros(N,N),Lc*D, zeros(N,N),zeros(N,N),zeros(N,M);
%    L^(-1)*B',zeros(M,N),zeros(M,N),zeros(M,N),-L^(-1)*R];

else


X=[-C_t^(-1)*Y_l,C_t^(-1),zeros(N,N),-C_t^(-1)*B; 
   alpha,beta,gamma,zeros(N,M);
   -eye(N),zeros(N),zeros(N),zeros(N,M); 
   L^(-1)*B',zeros(M,N),zeros(M,N),-L^(-1)*R];

[t,y]=ode45(@(t,y) ZIPloads(t,y,X,net_par,expo),tr,inval);
end





