function volt = ZIPloads(t,y,X,net_par,expo)

P_l   =  net_par.P_l;
Vref =  net_par.Vref;
C_t   =  net_par.C_t;
I_l  =  net_par.I_l;
M    =  net_par.M;
N    =  net_par.N;

I=[-C_t^(-1)*(I_l+P_l*y(1:N).^(expo-1));zeros(N,1);Vref;zeros(M,1)];

volt=X*y+I;

end