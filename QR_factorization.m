function [Q, R] = QR_factorization(A)

[rows ,col]=size(A);

%Matrix Q
Q(:,1) = A(:,1)/norm( A(:,1) );
for i=2:col
    q = A(:,i);
    for l=1:i-1
        q = q - ( A(:,i)' * Q(:,l) ) * Q(:,l) ;
    end
    
    Q(:,i) = q./norm(q) ;
end

%Matrix R
R=zeros(col,col);
R(1,1) = norm( A(:,1) );

for r=1:rows
    for c=2:col
        if(c>=r)
            R(r,c) = Q(:,r)' * A(:,c) ;
        end
    end
end
%max( max( abs(A-Q*R) ) )

end