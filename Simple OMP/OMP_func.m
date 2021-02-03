function [s_hat,t] = OMP_func(Fi,u,m)


%step 1
L = [];        %Index set of lt(lamda taf)
t = 0;         %iterations
Ft = [];
[N d] = size( Fi );
rm = u;

%step 1
for ii=1:m
    t=t+1;
    
    %step 2
    tmp = zeros(d,1);
    for i=1:d
        tmp(i,1) = abs( dot( rm,Fi(:,i) ) ) ; 
    end
    
    [M, lt] = max(tmp);
    
    %step 3
    %update index set
    L = [L; lt];
    
    %update matrix Ft
    Ft = [Ft Fi(:,lt)];
    

    %least squears problem
    %im not sure if Ft is non singular
    %o Ft einai full column rank opote to pinv einai koble
    xt = pinv(Ft) * u;     %to eipame me ton kurio liava

    %setp 5
    %Calculate the new approximation of data(they mean u) and the new residual
    am = Ft * xt;
    rm = u - am;
    
    %break when the residual's norm is smaller than a certain threshold
%     if( norm(rm) < 10^-8)
%         break;
%     end 
    
end

s_hat = zeros(d,1);
for i = 1:t %length( L(:,1) )
    s_hat( L(i,1),1 ) = xt(i);
end

end















