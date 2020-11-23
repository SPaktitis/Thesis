clc;
clear;

CSIT_com_nnoise =[];
CSIT_com_wnoise =[];
CSIT_real_nnoise =[];
T = 30:70;

for i=1:length(T)
    CSIT = JOMP_complex_noNoise(T(i));
    CSIT_com_nnoise = [CSIT_com_nnoise CSIT] ;
end

figure;
semilogy(T,CSIT_com_nnoise);
title("Complex channel without noise")
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");

for i=1:length(T)
    CSIT = JOMP_complex_wNoise(T(i));
    CSIT_com_wnoise = [CSIT_com_wnoise CSIT] ;
end

figure;
semilogy(T,CSIT_com_wnoise);
title("Complex channel with noise")
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");









