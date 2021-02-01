clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meleth kanaliwn apo statial optiki. Arxika create 2 kanalia kai valta na
% ftanoun sto dekti me diaforetikes gwnies, pai3e me tis times kai
% diapistwse pote o matrix H einai full rank kai pote oxi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lamda_c = .001; %signal's wavelength
d1 = 10^3; %distance between the antenna arrays
d2 = 2* 10^3;
Lr = 16;
Lt = 16;
nr = 32;
nt = 32;
Dt = Lt/nt;
Dr = Lr/nr;


%antenna gain
a1 = sqrt(.5) * ( randn + 1i*randn );
a2 = sqrt(.5) * ( randn + 1i*randn );

%Omegai = cos( phi[i] )
Omegat1 = cos( pi/2 );
Omegat2 = cos( .5111111111111111*pi );

Omegar1 = cos( pi/6 );
Omegar2 = cos( pi/5 );

%Creation of angular domain vectors et and er
et1 = Er(Omegat1,nt,Dt);
et2 = Er(Omegat2,nt,Dt);

er1 = Er(Omegar1,nr,Dr);
er2 = Er(Omegar2,nr,Dr);


%channel with arrays at the transmitter and the receiver with line of sight
%H = a * sqrt(nt*nr) * exp( -1i*2*pi*d/lamda_c ) * er * et' ;


% Channel matrix with 2 Tx far apart from each other and antenna array Rx
h1 = a1* sqrt(nr)* exp( -1i*2*pi*d1/lamda_c ) * er1;
h2 = a2* sqrt(nr)* exp( -1i*2*pi*d2/lamda_c ) * er2;

H = [h1, h2];
%rank(H)

%channel matrix with transmitter antenna array Tx and 2 antennas Rx far appart
h1 = a1 * sqrt(nt)* exp( -1i*2*pi*d1/lamda_c ) * et1;
h2 = a2 * sqrt(nt)* exp( -1i*2*pi*d1/lamda_c ) * et2;

H = [h1'; h2'];
rank(H)

%%%%%  Trying to reproduce the polarplot(figure 7.5)  %%%%%
Lr = 4;
nr = 8;
Dr = Lr/nr;
phi = 0:.01:2*pi

phi0 = pi/2;

er1 = Er(cos(phi0),nr,Dr);

vector=[];
for i = phi
   er2 =  Er(cos(i),nr,Dr) ;
   vector = [vector;  er1' * er2 ];
   
end


figure;
polarplot(phi,vector);
%ylim( [0 1] )
