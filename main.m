clc;
clear;

%variales
m=10;   %sparsity lvl
N=50;   
d=300;  

s=zeros(d,1); %arbitrary signal to recover
%----manualy creating the sparse signal---------
%s(1,1)=round(100*(rand(1,1)+1));
s(30:30:300,1)=round(100*(rand(m,1)+1));

Fi=normrnd(0,1/N,N,d);
u=Fi*s;

[s_hat,L,am,rm]=OMP(Fi,u,m);







