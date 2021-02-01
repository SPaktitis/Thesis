
clc;
clear;

d = 256;
m = [4 12 20 28 36];
N = 4:250;
figure;

for j = m
    F = [];
    percentage = [];
    for i = N
        count=0;
        for k=1:10^3
            s = zeros(d,1);
            indexes = randi([1 , d],j,1); %for easy check at sparsity indexes
            s(indexes,1) = randi([1, d],j,1);
            %Fi = normrnd(0,1/i,i,d);
            F = 1/sqrt(i) .* sign( randn(i,d) );
            u = F*s;
        
        
            [s_hat taph] = OMP_func(F,u,j);
            
            if( norm( s-s_hat )< 10^-3 )
                count = count +1;
            end
             
        end
        percentage = [percentage count/10];
    end
    
    plot(N,percentage,'-o');
    hold on;

end
hold off;





