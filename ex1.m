A = input("");
[L, U] = LU_decomposicao(A);
disp("L = ");
disp(L);
disp("U = ");
disp(U);

function [L,U] = LU_decomposicao(A)
n = size(A,1);
L = zeros(n);
U = eye(n);
    for k=1:n
        L(k,k) = A(k,k) - L(k,1:k-1) * U(1:k-1,k);
        for j=k+1:n
            U(k,j) = (A(k,j) - L(k,1:k-1) * U(1:k-1,j)) / L(k,k);
        end
        for i=k+1:n
            L(i,k) = (A(i,k) - L(i,1:k-1) * U(1:k-1,k)) / U(k,k);
        end
    end
end