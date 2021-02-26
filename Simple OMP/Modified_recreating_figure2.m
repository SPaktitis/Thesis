%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auto pou alazei edw einai pws to correct recovery tou signal ginetai
%   me vash th swsth anagnwrish twn mh mhdenikwn 8esewn


clc;
clear; 
%OMP tryouts with real variables
figure;
dim = 1;

sparsity = 0:5:80;
n = [52 100 148 196 244] ;

for jj = n
percent = [];
for kk = sparsity
count = 0;
taph  = [];  


for ii = 1:10^3
    
%variales
m = kk;          %kk;          %40;   %sparsity lvl
N = jj;          %75;         %200;   
d = 256;         %800;
s = zeros(d,dim);   %arbitrary signal to recover
%----manualy creating the sparse signal---------
indexes = randi([1 , d],m,1); %for easy check at sparsity indexes

s( indexes,1:dim ) = 1; %sqrt(.5).*( randn(m,dim) + 1i*randn(m,dim) );;%randi([1, 100],m,1);

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

for it = 1:kk
    %t=t+1;
    
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
    %im not sure if Ft is non singular
    %o Ft einai full column rank opote to pinv einai koble
    xt = pinv(Ft)*u;     %to eipame me ton kurio liava

    %setp 5
    %Calculate the new approximation of data(they mean u) and the new residual
    am = Ft * xt;
    rm = u - am;
    
%     %break when the residual's norm is smaller than a certain threshold
%     if( norm(rm) < 10^-8)
%         break;
%     end 
    
end

s_hat = zeros(d,dim);
for i=1:t
    s_hat(L(i,1),:) = xt(i,:);
end

taph = [taph; t];

%correctly recovered signal
test_vec = zeros(d,dim);
test_vec( L(:),: ) = 1;
if (  (s-test_vec) == 0 )
    count = count + 1;
end 


end

fprintf("success percentage: "+count/10+"   ");
fprintf("iterations percentage: "+sum(taph)/1000+"\n");

percent = [percent count/10];


end

plot(sparsity,percent,'-o');
hold on;
end
hold off;
legend('m=4','m=12','m=20','m=28','m=36');
ylabel("% recovered signals");
xlabel("Number of measurements (N)");
title("% of input signals recovered correctly (d=256)");
