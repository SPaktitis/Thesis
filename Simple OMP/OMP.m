clc;
clear; 
%OMP tryouts with real variables
figure;

M = [4 12 20 28 36];
n = 1:5:260 ;

for kk = M
percent = [];
for jj = n
count=0;
taph =[];  


for ii=1:10^3
    
%variales
m = kk;          %40;   %sparsity lvl
N = jj;    %75;         %200;   
d = 256;        %800;
s=zeros(d,1);   %arbitrary signal to recover
%----manualy creating the sparse signal---------
indexes = randi([1 , d],m,1); %for easy check at sparsity indexes

s(indexes,1) = 1;%randi([1, 100],m,1);

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
    for i=1:d
        tmp(i,1)=abs(dot(rm,Fi(:,i))); 
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
    if( norm(rm) < 10^-9)
        break;
    end 
    
end

s_hat = zeros(d,1);
for i=1:t
    s_hat(L(i,1),1) = xt(i);
end

taph = [taph; t];

if ( norm( s-s_hat )< 10^-5 )
    count = count + 1;
end 


end

fprintf("success percentage: "+count/10+"   ");
fprintf("iterations percentage: "+sum(taph)/1000+"\n");

percent = [percent count/10];


end

plot(n,percent);
hold on;
end
hold off;
legend('m=4','m=12','m=20','m=28','m=36');
ylabel("% recovered signals");
xlabel("Number of measurements (N)");
title("% of input signals recovered correctly (d=256)");