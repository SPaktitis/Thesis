function [s_hat,L,am,rm]=OMP(Fi,u,m)
%OMP tryouts with real variables


% s=zeros(d,1); %arbitrary signal to recover
% %----manualy creating the sparse signal---------
% %s(1,1)=round(100*(rand(1,1)+1));
% s(30:30:300,1)=round(100*(rand(m,1)+1));
% %-------------------------------------------
% 
% %Fo einai o pinakas pou dhmiourgw mesa sta loops...
% %making of the mesurement matrix Fi
% %with Normal Distribution
% Fi=normrnd(0,1/N,N,d);

%step 1
L=0;   %Index set of lt(lamda taf)
%u=Fi*s;
rm=u;    %arxikopoihsh residual 
t=0;    %iterations
while (1)
    t=t+1;
    %step 2
    tmp=zeros(d,1);
    for i=1:d
        tmp(i,1)=abs(dot(r,Fi(:,i))); 
    end
    %auto douleuei twra pou einai mono ena to megisto
    [M,lt]=max(tmp);
    
    %step 3
    %update index set
    if(t>1)
        L(t,1)=lt;
    else
        L=lt;
    end;
    %update matrix Fo
    Fo(:,t)=Fi(:,lt);

    %least squears problem
    %im not sure if Fo is non singular
    %o Fo einai full column rank opote to pinv einai koble
    xt=pinv(Fo)*u;     %ti eipame me ton kurio liava

    %setp 5
    %Calculate the new approximation of data(they mean u) and the new residual
    am=Fo*xt;
    rm=u-am;
    
    %reak when the residual's norm is less than 10^-5
    if(norm(rm)<10^-4)
        break;
    end    
end

s_hat=zeros(d,1);
for i=1:m
    s_hat(L(i,1),1)=xt(i);
end

















