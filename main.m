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
% axis([30 75 0.0001 1])
% title("NMSE of CSIT Versus T")
% xlabel("CSIT Estimation Overhead T");
% ylabel("NMSE of CSIT ");


%%%%%%%%% Different snr %%%%%%%%%%%%%%%%%
% clc;
% clear;
% 
% NMSE_diff_snr=[];
% snr = 15:2.5:35;
% for i = snr
%     CSIT = JOMP_diff_snr( i );
%     NMSE_diff_snr = [NMSE_diff_snr CSIT] ;
%     fprintf("SNR "+i+"dB, CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(snr,NMSE_diff_snr);
% title("NMSE of CSIT Versus SNR")
% axis([15 35 0.00001 1])
% xlabel("SNR(dB)");
% ylabel("NMSE of CSIT");


%%%%%%%%%%% Different common support sc %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
% 
% NMSE_diff_sc =[];
% sc = 2:10;
% for i=sc
%     CSIT = JOMP_diff_sc( i );
%     NMSE_diff_sc = [NMSE_diff_sc CSIT] ;
%     fprintf("Sc ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(sc,NMSE_diff_sc);
% axis([2 10 0.0000001 1])
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
% axis([14 20 0.0000001 1])
% title("NMSE of CSIT Versus S")
% xlabel("S (individual sparsity)");
% ylabel("NMSE of CSIT");

%%%%%%%%%%% Different number of users K %%%%%%%%%%%%%%%%%%%
% clc;
% clear;
%  
% NMSE_diff_K =[];
% K = 10:5:55;
% for i = K
%     CSIT = JOMP_diff_K( i );
%     NMSE_diff_K = [NMSE_diff_K CSIT] ;
%     fprintf("K ="+i+", CSIT= "+CSIT+"\n" );
% end
% 
% figure;
% semilogy(K,NMSE_diff_K);
% axis([10 55 0.0000001 1])
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
% axis([1 6 0.0001 1])
% title("NMSE of CSIT Versus N")
% xlabel("N (number of antennas at users)");
% ylabel("NMSE of CSIT");

%%%%%%%%%%%% Different number of antennas(M) at the base station %%%%%%%%%%%%%%%%%%%
clc;
clear;

NMSE_diff_M =[]; 
M = 60:20:220 ;
for i = M
    CSIT = JOMP_diff_M( i );
    NMSE_diff_M = [NMSE_diff_M CSIT] ;
    fprintf("M ="+i+", CSIT= "+CSIT+"\n" );
end

figure;
semilogy(M,NMSE_diff_M);
axis([60 220 0.0001 1])
title("NMSE of CSIT Versus N")
xlabel("N (number of antennas at users)");
ylabel("NMSE of CSIT");



