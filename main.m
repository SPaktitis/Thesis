%%%%%%%%%%%% Different snr %%%%%%%%%%%%%%%%%
clc;
clear;

NMSE_diff_snr=[];
snr = 15:35;
for i=1:length(snr)
    CSIT = JOMP_diff_snr(snr(i));
    NMSE_diff_snr = [NMSE_diff_snr CSIT] ;
end

figure;
semilogy(snr,NMSE_diff_snr);
title("NMSE of CSIT Versus SNR")
xlabel("SNR(dB)");
ylabel("NMSE of CSIT");

%%%%%%%%%%%% Different overhead T %%%%%%%%%%%%%%%%%
clc;
clear;


NMSE_diff_T =[];
T = 31:3:70;
for i=1:length(T)
    CSIT = JOMP_diff_T(T(i));
    NMSE_diff_T = [NMSE_diff_T CSIT] ;
end

figure;
semilogy(T,NMSE_diff_T);
title("NMSE of CSIT Versus T")
xlabel("CSIT Estimation Overhead T");
ylabel("NMSE of CSIT ");


%%%%%%%%%%%% Different common support sc %%%%%%%%%%%%%%%%%%%
clc;
clear;


NMSE_diff_sc =[];
sc = 2:10;
for i=1:length(sc)
    CSIT = JOMP_diff_sc(sc(i));
    NMSE_diff_sc = [NMSE_diff_sc CSIT] ;
end

figure;
semilogy(sc,NMSE_diff_sc);
title("NMSE of CSIT Versus Sc")
xlabel("Sc (common sparsity level parameter)");
ylabel("NMSE of CSIT");



