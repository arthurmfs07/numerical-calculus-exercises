% Display matrices

A = input("");

[L,U,P,P_matrices,M_matrices] = lup_decomp_with_permutation_matrices(A);
A_inverse = inv(U) * inv(L) * P;  % Inverse using LU decomposition

fprintf("inv(A) =\n");
disp(A_inverse);
fprintf("L =\n");
disp(L);
fprintf("U =\n");
disp(U);
fprintf("P =\n");
disp(P);
% Note that Pk and Mk have dimensions nxnxk, where k is the number of matrices
for k=1:size(A,1)-1
    fprintf("P%d =\n",k);
    disp(P_matrices(:,:,k)); % Pk contains all permutation matrices
end
for k=1:size(A,1)-1
    fprintf("M%d =\n",k);
    disp(M_matrices(:,:,k)); % Mk contains all elimination matrices
end


function [L,U,P,P_matrices,M_matrices] = lup_decomp_with_permutation_matrices(A)
% A: matriz nao-singular
% L, U: matrizes triang. inf. e sup., respectivamente
% P: matriz de permutacao
    n = size(A,1);
    P = eye(n); L = P; U = A;
    P_matrices = zeros(n,n,n-1); % Array to store permutation matrices
    M_matrices = zeros(n,n,n-1); % Array to store elimination matrices

    for k=1:n-1
        [~,p] = max(abs(U(k:n,k)));
        p = p+(k-1);
        P([k p],:) = P([p k],:);
        U([k p],k:n) = U([p k],k:n);
        L([k p],1:k-1) = L([p k],1:k-1);
        
        % Compute elimination matrix M_k
        M = eye(n);
        M(k+1:n, k) = -U(k+1:n, k) / U(k, k);
        M_matrices(:,:,k) = M;
        
        % Compute permutation matrix P_k
        P_k = eye(n);
        P_k([k, p], :) = P_k([p, k], :);
        P_matrices(:,:,k) = P_k;
        
        for i=k+1:n
            L(i,k) = U(i,k)/U(k,k);
            U(i,k:n) = U(i,k:n) - L(i,k)*U(k,k:n);
        end
    end
end

