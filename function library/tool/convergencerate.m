function B = convergencerate(A)

    [m,n] = size(A);
    B = zeros(m,2*n-1);
    B(:,1) = A(:,1);
    B(:,2) = A(:,2);
    for i = 4:2:2*n-1
        B(:,i) = A(:,i/2+1);
    end
    for i = 3:2:2*n-1
        for j = 2:m
            B(j,i) = log(B(j,i-1)/B(j-1,i-1))/log(B(j,1)/B(j-1,1));
        end
    end


end