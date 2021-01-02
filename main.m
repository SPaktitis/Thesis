%%%%%%%%%%% Different overhead T %%%%%%%%%%%%%%%%%
% clc;
% clear;
% 
% 
% NMSE_diff_T =[];
% T = 30:5:75;
% for i = T
%     CSIT = JOMP_diff_T( i );
%     NMSE_diff_T = [NMSE_diff_T CSIT] ;
%     fprintf("T= "+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(T,NMSE_diff_T);
% title("NMSE of CSIT Versus T")
% xlabel("CSIT Estimation Overhead T");
% ylabel("NMSE of CSIT ");


%%%%%%%%% Different snr %%%%%%%%%%%%%%%%%
clc;
clear;

NMSE_diff_snr=[];
snr = 15:2.5:35;
for i = snr
    CSIT = JOMP_diff_snr( i );
    NMSE_diff_snr = [NMSE_diff_snr CSIT] ;
    fprintf("SNR "+i+"dB, CSIT= "+CSIT+"\n" );
end

figure;
semilogy(snr,NMSE_diff_snr);
title("NMSE of CSIT Versus SNR")
xlabel("SNR(dB)");
ylabel("NMSE of CSIT");


%%%%%%%%%%% Different common support sc %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
%
% NMSE_diff_sc =[];
% sc = 2:1:10;
% for i=sc
%     CSIT = JOMP_diff_sc( i );
%     NMSE_diff_sc = [NMSE_diff_sc CSIT] ;
%     fprintf("Sc ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(sc,NMSE_diff_sc);
% title("NMSE of CSIT Versus Sc")
% xlabel("Sc (common sparsity level parameter)");
% ylabel("NMSE of CSIT");

%%%%%%%%%%% Different individual support S %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
%  
% NMSE_diff_s =[];
% s = 14:20;
% for i = s
%     CSIT = JOMP_diff_s( i );
%     NMSE_diff_s = [NMSE_diff_s CSIT] ;
%     fprintf("S ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(s,NMSE_diff_s);
% title("NMSE of CSIT Versus S")
% xlabel("S (individual sparsity)");
% ylabel("NMSE of CSIT");

%%%%%%%%%%% Different number of users K %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
%  
% NMSE_diff_K =[];
% K = 10:10:55;
% for i=1:length(K)
%     CSIT = JOMP_diff_K(K(i));
%     NMSE_diff_K = [NMSE_diff_K CSIT] ;
%     fprintf("K ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(K,NMSE_diff_K);
% title("NMSE of CSIT Versus K")
% xlabel("K (number of users)");
% ylabel("NMSE of CSIT");

%%%%%%%%%%%% Different number of antennas(N) at each user %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
% 
% NMSE_diff_N =[];
% n = 1:6;
% for i = n
%     CSIT = JOMP_diff_N( i );
%     NMSE_diff_N = [NMSE_diff_N CSIT] ;
%     fprintf("N ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(n,NMSE_diff_N);
% title("NMSE of CSIT Versus N")
% xlabel("N (number of antennas at users)");
% ylabel("NMSE of CSIT");



