clc;
clear;

M=160;              %transmit antennas
N=2;                %receive antennas
K=40;               %number of users
sc=4;               %common sparsity parameter
s=10;               %individual sparsity parameter
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
Xa = sqrt(P/M) .* (sign(2*rand(M,T)-1)) ;
X = At * Xa;

%Creation of the concatenated 
%Channel matrix Hw for K users
Hw=zeros(N*K,M);
Omegai = randi([1 M],K,s-sc);
Omegac=randi([1 M],sc,1);
for i=1:K   
   Hw(i*N-1:i*N , Omegai(i,:))  = randi([1, 100],2,length(Omegai(i,:)) ); 
   Hw(i*N-1:i*N , Omegac(:))    = randi([1, 100], 2,length(Omegac) );
end

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

%step2(Common support identification)
R = real( Y_hat );
Omegac_est = [];

for k = 1:sc
    
    paths = [];
    times = [];
    for j = 1:K %for each user calculate Omegai_est
        
        indexes=[];
        Ft = [];
        rm = R(:,j*N-1:j*N);
        
        %find the sc - |Omegac_est| columns we need
        it=0;
        while(1)
            it=it+1;    %iterations
            
            %Modified OMP to solve problem at A for 1 user            
            for l=1:M
                tmp1(l) = norm( X_hat(:,l)' *rm ) /norm(X_hat(:,l)) ; 
            end 
        
            [value , index] = max(tmp1);
            indexes = [indexes index];           
            Ft = [Ft X_hat(:,index)];       
            x2t = pinv(Ft) * R(:,j*N-1:j*N);
            at = Ft * x2t;
            rm = R(:,j*N-1:j*N) - at;
            
            if( norm(rm)<10^-6 )
                break
            end
        end     
                
        
        %====== B (Support pruning)
        l=[];
        for i=1:length(indexes)
            if(  norm(X_hat(:,indexes)' * R(:,j*N-1:j*N) , 'fro')^2 >= (eta1*N) ) 
                l = [l indexes(i)]; 
            end    
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
                   end %end if
               end %endfor paths matrix 
               
               % l(i) is not in paths matrix
               if( not(flag) )
                   paths = [paths l(i)];
                   times = [times 1];
               end    
            end %endfor l matrix            
        end %end if
    end %end for all users
    
    %======= C(Support Update)
    [value , index] = max(times);
    %bellow is a some code to deal with a situational problems
    %where the last support index is not retrieved correctly
    if( not(isempty(Omegac_est)) )
        t=1;
        while(1)
          if( paths(index) == Omegac_est(t) )
              times(index)   = 0;
              [value, index] = max(times);
              t=1;
          else
              t=t+1;
          end
          
          if( t>length(Omegac_est) )
              break;
          end
        end
       
    end    
    %update the indexes
    Omegac_est = [Omegac_est paths(index)];
    
    %======== D(Residual update)
    L = real( X_hat(:,Omegac_est) * pinv(X_hat(:,Omegac_est)) );
    
    for ii=1:K
       R(:,ii*N-1:ii*N ) = ( diag( ones( length(X_hat(:,1)) ,1) ) - L ) * Y_hat(:,ii*N-1:ii*N);
    end
       
end


%===================  STEP 3 ==================================
Omegai_est = {};
%R = real( Y_hat );
L =[];
Omega_vector =[];
for i=1:K %for all users
    Omega_vector = Omegac_est;
    t=0;    %iterations counter   
    while (1)
        t=t+1;
        %Alfa(Upport Update)       
        tmp=[];
        for j=1:M
            tmp(j) = norm( X_hat(:,j)' *R(:, i*N-1:i*N) ) /norm(X_hat(:,j));
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
        
        if( t >= (s - sc) )
            break;
        end
        
        
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
   
   H_est_hat(sort(vector) , i*N-1:i*N ) = pinv( X_hat(:,sort(vector) ) ) * Y_hat(: ,i*N-1:i*N) ;
   
   H_est(i*N-1:i*N,:) = Ar * H_est_hat(:,i*N-1:i*N)' * At' ;
   
end