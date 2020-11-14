clc;
clear;

M=160;              %transmit antennas
N=2;                %receive antennas
K=40;               %number of users
sc=4;               %common sparsity parameter
s=6;               %individual sparsity parameter
P=28;               %transmit SNR in dB
eta1=0.2;           %parameters used 
eta2=2;             %in JOMP alg.
Dt=1/2;             %antenna spacing
Dr=1/2;             %antennas spacing
Lt=round(M/2);      %Transmit antenna length 
Lr=round(N/2);      %Receive antenna length
T=50;               %number of pilot symbols



%Creation of angular basis matrix At and Ar
e=[];
At=[];
for k=1:M
    for a=1:M
        tmp1 = exp(-1i*2*pi*(a-1)*Dt*(k-1)/Lt);
        e = [e;    tmp1];
    end
    et = 1/sqrt(M) .* e;
    At = [At  et];
    e=[];
    et=[];
end

e=[];
Ar=[];
for k=1:N
    for a=1:N
        tmp1 = exp(-1i*2*pi*(a-1)*Dr*(k-1)/Lr);
        e = [e;    tmp1];
    end
    er = 1/sqrt(N) .* e;
    Ar = [Ar  er];
    e=[];
    er=[];
end
%OMP tryouts with real variables

%----manualy creating the sparse signal---------
Hw = zeros(N,M);
vector = randi([1 , M],1,sc+s); %for easy check at sparsity indexes

Hw(1:N,vector) = randi([1, 100],2,sc+s); % +1i*randi([1, 100],2,sc+s);

H = Ar * Hw * At ;

%-------------------------------------------
Xa = sqrt(P/M) .* (sign(2*rand(M,T)-1)) ;
X = At * Xa;

Y = Xa' * Hw' ;

%Calculate hat amounts
X_hat = sqrt(M/(P*T)) .* (X' * At);
H_hat = Hw' ;
Y_hat = X_hat * H_hat;


%step 1
L=[];               %Index set of lt(lamda taf)
u = Y;              %u is Y_hat TxN 50x2
rm = u;             %arxikopoihsh residual 
t=0;                %iterations
Ft=[];

while (1)
    t=t+1;
    
    %step 2
    %tmp = zeros(d,1);
    for i = 1 : M
        %tmp(i,1)=abs(dot(rm,Fi(:,i)));
        tmp1(i) = norm( Xa(i,:) * rm) / norm( Xa(i,:) );
    end
    
    [M,lt] = max(tmp1);
    
    %step 3
    %update index set
    L = [L; lt];
    
    %update matrix Ft
    Ft = [Ft Xa(lt,:)'];
    

    %least squears problem
    %im not sure if Ft is non singular
    %o Ft einai full column rank opote to pinv einai koble
    xt = pinv(Ft) * u;     %to eipame me ton kurio liava

    %setp 5
    %Calculate the new approximation of data(they mean u) and the new residual
    am = Ft * xt;
    rm = u - am;
    
    %break when the residual's norm is smaller than a certain threshold
    if( norm(rm) < 10^-6)
        break;
    end 
    
end

s_hat = zeros(d,1);
for i=1:t
    s_hat(L(i,1),1:2) = xt(i,:);
end
