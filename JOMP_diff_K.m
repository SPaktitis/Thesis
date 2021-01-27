function CSIT = JOMP_diff_K(Kappa)

M      =160;           %transmit antennas
N      =2;             %receive antennas
K      =Kappa;         %number of users
sc     =9;             %common sparsity parameter
s      =17;            %individual sparsity parameter
SNR_dB =28;            %transmit SNR in dB
eta1   =0.2;           %parameters used 
eta2   =2;             %in JOMP alg.
Dt     =1/2;           %antenna spacing
Dr     =1/2;           %antennas spacing
Lt     =round(M/2);    %Transmit antenna length 
Lr     =round(N/2);    %Receive antenna length
T      =45;            %number of pilot symbols

P = M * 10^(SNR_dB/10) ;    %quantity used to adjust the trasnmitt snr
NMSE=[];


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



for pkt_num=1:100
    

%Creation of the concatenated 
%Channel matrix Hw for K users
Hw      =  zeros(N*K,M);
Omegai  = {};
Omegac  = unique( randi([1 M],sc,1), 'sorted' );
for i=1:K
   %the 2 lines of code bellow are used to generate a sparsity value around the given
   %boundary with a small but completely specified variance
   temp = randi([1 M],1,randi([s-2 s],1,1) - randi([sc sc+2],1,1)) ;
   Omegai  = [Omegai;  {temp} ];
   
   Hw(i*N-(N-1):i*N , temp)  = sqrt(.5) * ( randn(N,length( temp )) +...
                                             1i *randn( N,length( temp ) ) );

   Hw(i*N-(N-1):i*N , Omegac(:))    = sqrt(.5) * ( randn(N,length(Omegac)) +...
                                           1i *randn( N,length(Omegac) ) );
end

%Creation of the concatenated channel matrix
%Hi for all K users
H = [];
for i=1:K
    H = [H; Ar * Hw(i*N-(N-1):i*N,:) * At' ];
end

%Pilot matrix X
Xa = sqrt(P/M) .* (sign(2*rand(M,T)-1)) ;
X = At * Xa;

%%%%%%%%%
Y = zeros(N*K,T);
Noise = sqrt(.5) .* ( randn( size(Y) ) + 1i *randn( size(Y) ) );
for i=1:K
   Y(i*N-(N-1):i*N,:) = H(i*N-(N-1):i*N,:) * X + Noise(i*N-(N-1):i*N,:) ; 
end 

%============== Beggining of the algorithm =========================


    %step1
    %Calculate hat amounts
    X_hat = sqrt(M/(P*T)) .* (X' *At);
    %H_hat = Hw' ;

    Y_hat=[];
    for j=1:K
        Y_hat(:,j*N-(N-1):j*N) = sqrt(M/(P*T)) .* ( Y(j*N-(N-1):j*N,:)' *Ar);
    end

    %N_hat = sqrt(M/(P*T)) .* N' *Ar;

    %step2(Common support identification)
    R = Y_hat ;
    Omegac_est = [];

    for k = 1:sc
    
        paths = [];
        times = [];
        for j = 1:K %for each user calculate Omegai_est
        
            indexes=[];
            Ft = [];
            rm = R(:,j*N-(N-1):j*N);
        
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
                x2t = pinv(Ft) * R(:,j*N-(N-1):j*N);
                at = Ft * x2t;
                rm = R(:,j*N-(N-1):j*N) - at;
            
                if( norm(rm)<10^-6 )
                    break
                end
            end     
                
        
            %====== B (Support pruning)
            l=[];
            for i=1:length(indexes)
                if(  norm(X_hat(:,indexes)' * R(:,j*N-(N-1):j*N) , 'fro')^2 >= (eta1*N) ) 
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
        L = real( X_hat(:,Omegac_est) * pinv(X_hat(:,Omegac_est)) ) ;
    
        for ii=1:K
           R(:,ii*N-(N-1):ii*N ) = ( diag( ones( length(X_hat(:,1)) ,1) ) - L ) * Y_hat(:,ii*N-(N-1):ii*N);
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
                tmp(j) = norm( X_hat(:,j)' *R(:, i*N-(N-1):i*N) ) /norm(X_hat(:,j));
            end
            
            [value, index] = max(tmp);
            
            Omega_vector = [Omega_vector index];
        
        
            %B(Residual Update)
            L = real( X_hat(:,Omega_vector) * pinv( X_hat(:,Omega_vector) ) );
        
            R(:, i*N-(N-1):i*N) = (diag(ones(length(X_hat(:,1)) ,1)) - L ) *Y_hat(:, i*N-(N-1):i*N);
    
    
            %terminating conditions
            if( norm(R(:, i*N-1:i*N), 'fro')^2 <= (eta2 * N *M)/P )
                break;
            end
        
            if( t >= (s - sc) )
                break;
            end 
        end %end while  
    
        %update Omega-_est matrix
        Omegai_est = [Omegai_est; {Omega_vector}];
        
    
    end %end for

    
    %============== STEP4 ===================
    H_est = zeros(N*K,M);
    H_est_hat = zeros(M,N*K);
    for i=1:K
        vector = [];
   
        vector = cell2mat(Omegai_est(i,1));
   
        H_est_hat(sort(vector) , i*N-(N-1):i*N ) = pinv( X_hat(:,sort(vector) ) ) * Y_hat(: ,i*N-(N-1):i*N) ;
   
        H_est(i*N-(N-1):i*N,:) = Ar * H_est_hat(:,i*N-(N-1):i*N)' * At' ;
   
    end


    %=========== NMSE
    NMSE= norm( H - H_est, 'fro' ).^2 / norm( H, 'fro' ).^2;
end
    CSIT = sum(NMSE)/pkt_num;

end