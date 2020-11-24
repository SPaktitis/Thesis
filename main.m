clc;
clear;

NMSE_com =[];
NMSE_com_n =[];
NMSE_real =[];
T = 31:3:70;


for i=1:length(T)
    CSIT = JOMP_real(T(i));
    CSIT_real = [NMSE_real CSIT] ;
end

figure;
semilogy(T,NMSE_real);
title("Real channel without noise")
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");

%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(T)
    CSIT = JOMP_com(T(i));
    CSIT_com = [NMSE_com CSIT] ;
end

figure;
semilogy(T,NMSE_com);
title("Complex channel without noise")
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(T)
    CSIT = JOMP_c_n(T(i));
    NMSE_com_n = [NMSE_com_n CSIT] ;
end

figure;
semilogy(T,NMSE_com_n);
title("Complex channel with noise")
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");









