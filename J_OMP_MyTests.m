clc;
clear;

M=160;              %transmit antennas
N=2;                %receive antennas
K=40;               %number of users
sc=9;               %common sparsity parameter
s=17;               %individual sparsity parameter
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

%Pilot matrix X
Xa = sqrt(P/M) .* (sign(2*rand(M,T)-1));
X = At * Xa;



%Creation of the concatenated 
%Channel matrix Hw for K users
Hw=zeros(N*K,M);
Omegai=[];
Omegac=randi([1 M],sc,1);
for i=1:K
   Omegai(i,:) = randi([1 M],s,1);
   Hw(i*N-1:i*N , Omegac(:))    = 1; 
   Hw(i*N-1:i*N , Omegai(i,:))  = 1;  
end

%altarnative method for Hw??
% Hw=[];
% for i=1:K
%     tmp=zeros(N,M);
%     Omegai(i,:) = randi([1 M],s-sc,1);
%     tmp(:,Omegai(i,:)) = 1;
%     tmp(:,Omegac) = 1;
%     Hw = [Hw; tmp];
% end


%Creation of the concatenated channel matrix
%Hi for all K users
H = [];
for i=1:K
    H = [H; Ar * Hw(i*N-1:i*N,:) *At ];
end
%===========================================
%Beggining of the algorithm

%step1
%Calculate hat amounts
X_hat = sqrt(M/(P*T)) .* (X' *At);
H_hat = Hw' ;
Y_hat = X_hat * H_hat;
%N_hat = sqrt(M/(P*T)) .* N' *Ar;


R = abs(Y_hat); %periexei migadika skoupidia
Omegac_est = [];


%==================== STEP 2 ====================================
for k = 1:sc
    
    %h logikh auth pisteuw douleuei giati o pinakas X kai X_hat einai full rank

    Omegai_est={};
    paths = [];
    times = [];
    for j = 1:K %for each user calculate Omegai_est
        
        
        
        %===============Alpha
        tmp=[];
        for i=1:M
            tmp(i,1) = norm((X_hat(:,i)' *R(:, j*N-1:j*N) ), 'fro');  %abs(dot(rm,Fi(:,i))); 
        end
        
        value=0;
        index=0;
        indexes=[];
        %find the sc - |Omegac_est| columns we need
        for i=1 : length(Omegai(j,:))-length(Omegac)
            [value,index]=max(tmp);
            indexes = [indexes index];
            tmp(index) = 0;   %gia na mhn vriskei to idio megisto olh thn wra
        end        
        %h timh ths posothtas sto "A" vhma
        %norm( (X_hat(:,sort(Omega_i_est))' *R), 'fro') 
    
        %====== B (Support pruning)
        l=[];
        for i=1:length(indexes)
            if(  norm(X_hat(:,indexes)' *R , 'fro')^2 >= (eta1*N) ) 
                l = [l indexes(i)]; 
            end    
        end
        
        %Update Omegai_est matrix
        if ( isempty(Omegai_est) )
            Omegai_est = {l};
        else
            Omegai_est = [Omegai_est; {l}];
        end
        
        %the following if/else is implementing the step C(Support Update)
        %update paths and times matrices
        if( isempty(paths) )
            paths = l;
            times = ones( 1,length(paths) );
        else
            for i=1:length(l) %for each element in vector l
               flag = 0;
               for ii = 1:length(paths)
                   if ( paths(1,ii) == l(i) ) %if the path allready exists
                       times(1,ii) = times(1,ii)+1;
                       flag =1;
                   end %endif
               end %endfor paths matrix
               
               % not(flag) := element l(i) is not in paths matrix
               if( not(flag) )
                   paths = [paths l(i)];
                   times = [times 1];
               end    
            end %endfor l matrix            
        end 
    end
    %======= C(Support Update)
    [value , index] = max(times);
    Omegac_est = [Omegac_est paths(index)];
    
    L =[];
    %======== D(Residual update)
    L = X_hat(:,Omegac_est) *pinv(X_hat(:,Omegac_est));
    
    for ii=1:K
       R( :, ii*N-1:ii*N ) = ( diag( ones( length(X_hat(:,1)) ,1) ) - L ) * Y_hat(:,ii*N-1:ii*N);
    end
    
    
end

%===================  STEP 3 ==================================
Omegai_est = {};
R = Y_hat;
L =[];
Omega_vector =[];
for i=1:K %for all users
    t=1;    %iterations counter
    while (1)
        %Alfa(Upport Update)
        Omega_vector = Omegac_est;
        
        tmp=[];
        for j=1:M
            tmp(j) = norm( (X_hat(:,j)' * R(:, i*N-1:i*N)), 'fro' );
        end
        [value, index] = max(tmp);
    
        Omega_vector = [Omega_vector index];
        
        
        %B(Residual Update)
        L = X_hat(:,Omega_vector) * pinv( X_hat(:,Omega_vector) );
        
        R(:, i*N-1:i*N) = (diag(ones(length(X_hat(:,1)) ,1)) - L ) *Y_hat(:, i*N-1:i*N);
    
    
        %terminating conditions
        if( norm(R(:, i*N-1:i*N), 'fro')^2 <= (eta2 * N *M)/P )
            break;
        end
        
        if( t >= (length(Omegai(i,:)) - length(Omegac)) )
            break;
        end
        
        t=t+1;
    end %end while  
    
    %update Omega-_est matrix
    if( isempty(Omegai_est) )
        Omegai_est = {Omega_vector};
    else
        Omegai_est = [Omegai_est; {Omega_vector}];
    end
    
end %end for
    
%============== STEP4 ===================

H_est = zeros(N*K,M);
H_est_hat = zeros(M,N*K);
for i=1:K
   vector = [];
   
   vector = cell2mat(Omegai_est(i,1));
   
   H_est_hat(vector , i*N-1:i*N ) = abs (pinv( X_hat(:,vector) ) * Y_hat(: ,i*N-1:i*N) );
   
   H_est(i*N-1:i*N,:) = Ar * H_est_hat(:,i*N-1:i*N)' * At' ;
   
end






