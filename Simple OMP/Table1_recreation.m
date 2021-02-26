%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auto pou alazei edw einai pws to correct recovery tou signal ginetai
%   me vash th swsth anagnwrish twn mh mhdenikwn 8esewn


clc;
clear; 

dim = 1;
sparsity = [5 10 15];  %[4 8 12 16 20];    
n = 0:5:260 ;

%code to help recreate table 1 from the OMP paper
vector = zeros( size( sparsity ) );
t=1;    %pointer
d = 1024;         %800;

for kk = sparsity
percent = [];
jj = 0;
while( jj<= 260 )
count = 0;
taph  = [];  


for ii = 1:10^3
    
%variales
m = kk;       %kk;          %40;   %sparsity lvl
N = jj;          %75;         %200;   

s = zeros(d,dim);   %arbitrary signal to recover
%----manualy creating the sparse signal---------
indexes = randi([1 , d],m,1); %for easy check at sparsity indexes

s( indexes,1:dim ) = 1; %sqrt(.5).*( randn(m,dim) + 1i*randn(m,dim) );;%randi([1, 100],m,1);

%-------------------------------------------

%Ft einai o pinakas pou dhmiourgw mesa sta loops...
%making of the mesurement matrix Fi
%with Normal Distribution
%Fi = normrnd(0,1/N,N,d);

% measurements matrix from (kind of)Bernouli distribution
Fi = sign(2*rand(N,d)-1);


%step 1
L  = [];        %Index set of lt(lamda taf)
u  = Fi*s;
rm = u;         %arxikopoihsh residual 
Ft = [];

for it = 1:kk
    
    %step 2
    tmp = zeros(d,1);
    for i = 1:d
        tmp(i,1)=   abs( dot( rm,Fi(:,i) ) ); %norm( Fi(:,i)' *rm ) /norm( Fi(:,i) );   %norm( X_hat(:,l)' *rm ) /norm(X_hat(:,l))
    end
    
    [M,lt] = max(tmp);
    
    %step 3
    %update index set
    L = [L; lt];
    
    %update matrix Ft
    Ft = [Ft Fi(:,lt)];
    

    %least squears problem
    %Ft is full column rank so inversion is a valid operation
    xt = pinv(Ft)*u;     

    %setp 5
    %Calculate the new approximation of data(data=u) and the new residual
    am = Ft * xt;
    rm = u - am;
    
end

s_hat = zeros(d,dim);
for i=1:t
    s_hat(L(i,1),:) = xt(i,:);
end

%correctly recovered signal
test_vec = zeros(d,dim);
test_vec( L(:),: ) = 1;
if (  (s-test_vec) == 0 )
    count = count + 1;
end 


end

fprintf("success percentage: "+count/10+"\n");


percent = [percent count/10];

%code to help recreate table 1 from the OMP paper
if( percent(end) >= 95 )
    %store the first value for which the desired rate is achieved
    vector(t) = jj;
    t=t+1;
    %break the loop an continue to the next value of N
    jj = 261;
end

%increament of the iteration counter
if( percent(end) < 90 )
    jj = jj + 5;
else
    jj = jj + 1;
end

end


end