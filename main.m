clc;
clear;

CSIT_errors =[];
T = 30:70;

for i=1:length(T)
    CSIT = JOMP_complexHw_noNoise(T(i));
    CSIT_errors = [CSIT_errors CSIT] ;

end

figure;
semilogy(T,CSIT_errors);
xlabel("Overhead T");
ylabel("NMSE of CSIT estimation");






