clc;
clear; 
%OMP tryouts with real variables
figure;
dim = 1;

sparsity = 1:10:400;
n = [52 100 148 196 244];

for jj = n
percent = [];

for kk = sparsity  
count = 0;
taph  = [];  


for ii = 1:500
    
%variales
m = kk;       %kk;          %40;   %sparsity lvl
N = jj;          %75;         %200;   
d = 256;         %800;
s = zeros(d,dim);   %arbitrary signal to recover
%----manualy creating the sparse signal---------
indexes = randi([1 , d],m,1); %for easy check at sparsity indexes

s( indexes,1:dim ) = sqrt(.5).*( randn(m,dim) + 1i*randn(m,dim) );%1;%randi([1, 100],m,1);

%-------------------------------------------

%Ft einai o pinakas pou dhmiourgw mesa sta loops...
%making of the mesurement matrix Fi
%with Normal Distribution
Fi = normrnd(0,1/N,N,d);


%step 1
L  = [];        %Index set of lt(lamda taf)
u  = Fi*s;
rm = u;         %arxikopoihsh residual 
t  = 0;         %iterations
Ft = [];

while (1)
    t=t+1;
    
    %step 2
    tmp = zeros(d,1);
    for i = 1:d
        tmp(i,1)=  abs( dot( rm,Fi(:,i) ) );    %norm( Fi(:,i)' *rm ) /norm( Fi(:,i) );   %norm( X_hat(:,l)' *rm ) /norm(X_hat(:,l))
    end
    
    [M,lt] = max(tmp);
    
    %step 3
    %update index set
    L = [L; lt];
    
    %update matrix Ft
    Ft = [Ft Fi(:,lt)];
    

    %least squears problem
    %im not sure if Ft is non singular
    %o Ft einai full column rank opote to pinv einai koble
    xt = pinv(Ft)*u;     %to eipame me ton kurio liava

    %setp 5
    %Calculate the new approximation of data(they mean u) and the new residual
    am = Ft * xt;
    rm = u - am;
    
    %break when the residual's norm is smaller than a certain threshold
    if( norm(rm) < 10^-8)
        break;
    end 
    
end

s_hat = zeros(d,dim);
for i=1:t
    s_hat(L(i,1),:) = xt(i,:);
end

%to monitor the iterations
taph = [taph; t];

if ( norm( s-s_hat )^2 < 10^-6 )
    count = count + 1;
end 


end

fprintf("success percentage: "+count/5+"   ");
fprintf("iterations percentage: "+sum(taph)/500+"\n");

percent = [percent count/10];


end

plot(sparsity,percent,'-o');
hold on;
end
hold off;
legend('N=52','N=100','N=148','N=196','N=244');
ylabel("% recovered signals");
xlabel(" Sparsity level (m) ");
title("% of input signals recovered correctly (d=256)");
