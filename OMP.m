clc;
clear;
%OMP tryouts with real variables

%variales
m=40;   %sparsity lvl
N=200;   
d=800;
s=zeros(d,1); %arbitrary signal to recover
%----manualy creating the sparse signal---------
vector = randi([1 , d],m,1); %for easy check at sparsity indexes

s(vector,1) = randi([1, 100],m,1);

%-------------------------------------------

%Ft einai o pinakas pou dhmiourgw mesa sta loops...
%making of the mesurement matrix Fi
%with Normal Distribution
Fi = normrnd(0,1/N,N,d);


%step 1
L=[];        %Index set of lt(lamda taf)
u = Fi*s;
rm=u;       %arxikopoihsh residual 
t=0;        %iterations
Ft=[];

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













